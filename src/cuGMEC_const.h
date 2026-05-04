/*
 * cuGMEC
 *
 * Copyright (C) 2025 Shiyang Liu
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#pragma once

#include <unistd.h>
#include <cstdint>
#include <type_traits>

static uint64_t getHostHash(const char* string) {
    // Based on DJB2a, result = result * 33 ^ char
    uint64_t result = 5381;
    for (int c = 0; string[c] != '\0'; c++) {
        result = ((result << 5) + result) ^ string[c];
    }
    return result;
}

static void getHostName(char* hostname, int maxlen) {
    gethostname(hostname, maxlen);
    for (int i = 0; i < maxlen; i++) {
        if (hostname[i] == '.') {
            hostname[i] = '\0';
            return;
        }
    }
}

/*--------------------------------Constant--------------------------------*/

const double QE = 1.6021766208e-19;
const double MP = 1.672621637e-27;
const double PI = 3.1415926535897932;
const double MU0 = 4.0 * PI * 1.0e-7;
const double KEV = 1000.0 * QE;
enum picType { Ion, Alpha, Beam };
enum disType { Maxwell, Slowing0, Slowing1, Slowing2, Slowing3 };
enum spaceType { spaceReal, spaceUniform };
enum velocityType { velocityReal, velocityUniform };
enum matrixType { Laplacian, Resistive, Perp2Phi, Perp2dNe, Perp2dTe, Perp2dPi, Perp2dPa, Perp2dPb };
using trueType = std::integral_constant<bool, true>;
using falseType = std::integral_constant<bool, false>;
template <typename... Ts>
inline constexpr bool allTrue = (std::is_same_v<Ts, trueType> && ...);
template <typename T>
inline constexpr bool toBool = std::is_same_v<T, trueType>;