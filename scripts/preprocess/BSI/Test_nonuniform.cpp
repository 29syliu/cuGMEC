/* This is a C++ MEX file for MATLAB.
C++ B-spline interpolation interface
https://github.com/12ff54e/BSplineInterpolation
*/

#include <iostream>
#include <array>
#include <tuple>
#include <variant>
#include <Interpolation.hpp>
#include "mex.hpp"
#include "mexAdapter.hpp"

template<typename T> struct make_inverse_index_sequence_impl;
template<std::size_t... Is> struct make_inverse_index_sequence_impl<std::index_sequence<Is...>>
{
    using type = std::index_sequence<sizeof...(Is) - Is - 1 ...>;
};
template<std::size_t N> using make_inverse_index_sequence = typename make_inverse_index_sequence_impl<std::make_index_sequence<N>>::type;


using namespace matlab::engine;
using namespace  matlab::data;

decltype(auto) get_matlab_arr(auto&& x, auto i, auto... is)
{
    if constexpr (sizeof...(is) == 0)
        return x[i];
    else
        return get_matlab_arr(x[i], is...);
}

template<std::size_t I = 0, typename F, std::size_t N>
void range_invoke_row_major(F&& f, const std::array<std::size_t, N>& range, auto... is)
{
    if constexpr (I == N)
        f(is...);
    else
    {
        for (std::size_t i = 0; i < range[I]; i++)
            range_invoke_row_major<I + 1>(std::forward<F>(f), range, is..., i);
    }
}

template<std::size_t I = 0, typename F, std::size_t N>
void range_invoke_col_major(F&& f, const std::array<std::size_t, N>& range, auto... is)
{
    if constexpr (I == N)
        f(is...);
    else
    {
        //std::cout << "range " << I << " dim " << range[N - 1 - I] << std::endl;
        for (std::size_t i = 0; i < range[N - 1 - I]; i++)
            range_invoke_col_major<I + 1>(std::forward<F>(f), range, i, is...);
    }
}

decltype(auto) reverse_invoke(auto&& f, auto&&... is)
{
    const auto tup = std::tuple{ is... };
    return[&] <std::size_t... Is>(std::index_sequence<Is...>)->decltype(auto) {
        return f(std::get<Is>(tup)...);
    }(make_inverse_index_sequence<sizeof...(is)>{});
}

using non_uniform_range_t = std::pair<TypedIterator<const double>, TypedIterator<const double>>;
using uniform_range_t = std::pair<double, double>;
using range_variant_t = std::variant<uniform_range_t, non_uniform_range_t>;

range_variant_t get_range(const TypedArray<double> range_i)
{
    //std::cout << "in function get range" << std::endl;
    if (std::max(range_i.getDimensions()[0], range_i.getDimensions()[1]) == 2) {
        std::cout << "uniform range" << std::endl;
        return uniform_range_t(double{ range_i[0] }, double{ range_i[1] });
    }
    else {
        std::cout << "non-uniform range" << std::endl;
        return non_uniform_range_t(range_i.begin(), range_i.end());
    }
}


class MexFunction : public matlab::mex::Function {

    ArrayFactory factory;
    std::shared_ptr<MATLABEngine> matlabPtr = getEngine();

public:
    MexFunction()
    {
    }


    //template<std::size_t Dim>
    //void interpolation(const uint64_t order, const TypedArray<bool>& is_periodic, const TypedArray<double>& range, const TypedArray<double>& mesh_in, const TypedArray<double>& coor_in, const std::vector<TypedArray<std::size_t>>& derivative_vector, matlab::mex::ArgumentList& outputs)
    //{
    //    auto interpolation_function = interpolation_initial<Dim>(order, is_periodic, range, mesh_in);
    //    TypedArray<double> result = factory.createArray<double>({ coor_in.getDimensions()[0] });
    //    for (std::size_t i = 0; i < derivative_vector.size(); i++) {
    //        interpolation_nd<Dim>(interpolation_function, coor_in, derivative_vector[i], result);
    //        outputs[i] = result;
    //    }
    //}

    template<std::size_t Dim>
    void interpolation(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs)
    {
        if (inputs.size() == 6)
        {
            checkArguments(outputs, inputs);
            // const std::size_t task_num = inputs.size() - 5;
            const uint64_t order = inputs[0][0];
            const TypedArray<bool> is_periodic = std::move(inputs[1]);
            //const TypedArray<double> range = std::move(inputs[2]);
            const CellArray range = std::move(inputs[2]);
            //const TypedArray<double> range1 = range[0][0];
            const TypedArray<double> mesh_in = std::move(inputs[3]);
           // interpolation_initial<Dim>(order, is_periodic, range, mesh_in);

            const std::array<std::size_t, Dim> Dims = std::array<std::size_t, Dim>{mesh_in.getDimensions()[0]};
            auto f_nd = intp::Mesh<double, Dim>{ Dims[0] + std::size_t{ is_periodic[0] } };
            //for (std::size_t i = 0; i < Dims[0]; i++)
            //    f_nd[i] = mesh_in[i];
            range_invoke_col_major([&](auto... is)
                {
                    reverse_invoke(f_nd, is...) = get_matlab_arr(mesh_in, is...);
                }, Dims);

            using range_t = std::pair<TypedIterator<const double>, TypedIterator<const double>>;

            //range_t range_i = std::make_pair(range1.begin(), range1.end());
            //auto interp_func = intp::InterpolationFunction<double, Dim>{
            //    order, { is_periodic[0] }, f_nd, range_i };

            //auto range_unif = non_uniform_range_t(range1.begin(), range1.end());
            //auto range_var = range_variant_t{ range_unif };
           // auto range_x = get_range(range1);

            auto interp_func = std::visit([&](const auto& range_1)
                {
                    return intp::InterpolationFunction<double, Dim>{
                        order, { is_periodic[0] }, f_nd, range_1};
                }, get_range(range[0][0]));

            //checkArguments_interp(outputs, inputs);
            const TypedArray<double> coor_in = std::move(inputs[4]);
            const TypedArray<std::size_t> derivative_in = std::move(inputs[5]);
            TypedArray<double> result = factory.createArray<double>({ coor_in.getDimensions()[0] });
            //interpolation_nd<Dim>(coor_in, derivative_in, result);

            for (std::size_t i = 0; i < coor_in.getDimensions()[0]; i++)
            {
                result[i] = [&]<std::size_t... Is>(std::index_sequence<Is...>) {
                    return  interp_func.derivative({ double{ coor_in[i][Is]}... }, { derivative_in[Is]... });
                }(make_inverse_index_sequence<Dim>{});
            }

            outputs[0] = std::move(result);
        }
    }

    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) override {
        // inputs : ( order, is periodic, range, array, array interpolation, derivative )
        // outputs : ( result )

        if (inputs.size() < 2)
            matlabPtr->feval(u"error", 0,
                std::vector<Array>({ factory.createScalar("Two input at least required") }));
        std::size_t dim = inputs[1].getDimensions()[0];
        //std::cout << "dim " << dim << std::endl;

        switch (dim)
        {
        case 1:
            interpolation<1>(outputs, inputs);
            break;
        default:
            std::cout << "unsupport dim, you need add dim " << dim << " in bspline.cpp." << std::endl;
            throw;
        }
    }

    void checkArguments(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {

        //if (inputs.size() < 6) {
        //    matlabPtr->feval(u"error", 0,
        //        std::vector<Array>({ factory.createScalar("Six input at least required") }));
        //}

        if (inputs[0].getType() != ArrayType::UINT64 || inputs[1].getType() != ArrayType::LOGICAL || inputs[2].getType() != ArrayType::CELL || inputs[3].getType() != ArrayType::DOUBLE) {
            matlabPtr->feval(u"error", 0,
                std::vector<Array>({ factory.createScalar("Input type error") }));
        }

        const TypedArray<bool> is_periodic = inputs[1];
       // const TypedArray<double> range = inputs[2];
        const TypedArray<double> mesh_in = inputs[3];
        std::size_t dim = 0;
        for (std::size_t i = 0; i < mesh_in.getDimensions().size(); i++)
            if (mesh_in.getDimensions()[i] > 1) dim++;

        //if (range.getDimensions().size() != 2 || range.getDimensions()[0] != dim || range.getDimensions()[1] != 2) {
        //    matlabPtr->feval(u"error", 0,
        //        std::vector<Array>({ factory.createScalar("Range require Dim * 2 array") }));
        //}

        if (mesh_in.getDimensions()[0] == 1)
        {
            matlabPtr->feval(u"error", 0,
                std::vector<Array>({ factory.createScalar("Mesh should be col vector for 1 dimension") }));
        }

        if (is_periodic.getDimensions()[0] != dim)
        {
            matlabPtr->feval(u"error", 0,
                std::vector<Array>({ factory.createScalar("is_periodic require length dim array") }));
        }
    }

    void checkArguments_interp(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
        if (inputs[0].getType() != ArrayType::DOUBLE || inputs[1].getType() != ArrayType::UINT64) {
            matlabPtr->feval(u"error", 0,
                std::vector<Array>({ factory.createScalar("Input type error") }));
        }
        const TypedArray<double> coor_in = inputs[0];
        const TypedArray<std::size_t> derivative_in = inputs[1];
        std::size_t dim = coor_in.getDimensions()[1];
        if (coor_in.getDimensions().size() != 2)
        {
            matlabPtr->feval(u"error", 0,
                std::vector<Array>({ factory.createScalar("Interpolate coordinate require N * dim array") }));
        }
        if (derivative_in.getDimensions()[0] != dim)
        {
            matlabPtr->feval(u"error", 0,
                std::vector<Array>({ factory.createScalar("Interp coordinate dimension is not equal to derivative dimension") }));
        }
        if (derivative_in.getType() != ArrayType::UINT64)
        {
            matlabPtr->feval(u"error", 0,
                std::vector<Array>({ factory.createScalar("derivative type error") }));
        }
    }
    // inputs[4].getType() != ArrayType::DOUBLE

};
