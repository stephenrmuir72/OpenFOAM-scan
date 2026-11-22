/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Keysight Technologies
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Description
    Minimal stand-in for compilers without <execution> (eg, gcc-7):
    sequenced policy and detection mechanisms. [issue #3661]

    It is not entirely legitimate to inject these into std/execution
    namespaces, but since <execution> was unavailable it doesn't really
    conflict either.

    OpenFOAM issue: #3661

\*---------------------------------------------------------------------------*/

#ifndef Foam_ListPolicy_legacy_H
#define Foam_ListPolicy_legacy_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Normally included by ListPolicy.H directly
#include "ListPolicy.H"

#if __has_include(<execution>)
    #include <execution>
#else
    #define Foam_execution_fallback
#endif

#ifndef __cpp_lib_execution
#ifndef Foam_execution_fallback
#define Foam_execution_fallback
#endif
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Fallbacks for older compilers, incomplete installations

#ifdef Foam_execution_fallback

#include <iterator>
#include <utility>

//
// Minimal stand-in for defining/detecting execution policies (#3661)
//
namespace std
{
namespace execution
{
    struct sequenced_policy {};
    inline constexpr sequenced_policy seq{};
}

template<class T>
struct is_execution_policy : false_type {};

template<>
struct is_execution_policy<execution::sequenced_policy> : std::true_type {};

template<class T>
constexpr bool is_execution_policy_v
    = std::is_execution_policy<std::decay_t<T>>::value;

} // End namespace std

#endif  /* [Foam_execution_fallback] */


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Fallbacks for older compilers, incomplete installations

#ifdef Foam_execution_fallback

#include <algorithm>

// Minimal stand-in for stripping out execution policies from common
// std algorithms. Mostly a "wack-a-mole" to cover whatever is currently
// being used in the code base.
//
// NB: some algorithms (eg, std::reduce) may be missing or buggy even in
// their plain form in older installations so they are bypassed entirely
// for these shims.

namespace std
{

// <algorithm>

template<class Policy, class... Args>
inline auto copy(Policy&&, Args&&... args) -> std::enable_if_t
<
    std::is_execution_policy_v<Policy>,
    decltype(std::copy(std::forward<Args>(args)...))
>
{
    return std::copy(std::forward<Args>(args)...);
}


template<class Policy, class... Args>
inline auto fill(Policy&&, Args&&... args) -> std::enable_if_t
<
    std::is_execution_policy_v<Policy>,
    void
>
{
    std::fill(std::forward<Args>(args)...);
}


template<class Policy, class... Args>
inline auto fill_n(Policy&&, Args&&... args) -> std::enable_if_t
<
    std::is_execution_policy_v<Policy>,
    decltype(std::fill_n(std::forward<Args>(args)...))
>
{
    return std::fill_n(std::forward<Args>(args)...);
}


template<class Policy, class... Args>
inline auto for_each(Policy&&, Args&&... args) -> std::enable_if_t
<
    std::is_execution_policy_v<Policy>,
    void
>
{
    std::for_each(std::forward<Args>(args)...);
}


template<class Policy, class... Args>
inline auto transform(Policy&&, Args&&... args) -> std::enable_if_t
<
    std::is_execution_policy_v<Policy>,
    decltype(std::transform(std::forward<Args>(args)...))
>
{
    return std::transform(std::forward<Args>(args)...);
}


template<class Policy, class... Args>
inline auto sort(Policy&&, Args&&... args) -> std::enable_if_t
<
    std::is_execution_policy_v<Policy>,
    void
>
{
    std::sort(std::forward<Args>(args)...);
}


template<class Policy, class... Args>
inline auto stable_sort(Policy&&, Args&&... args) -> std::enable_if_t
<
    std::is_execution_policy_v<Policy>,
    void
>
{
    std::stable_sort(std::forward<Args>(args)...);
}

// <numeric>
//
// Many of these may be missing entirely (even the serial version)
// so supply suitable replacements here, but wrapped with the execution
// policy to avoid any potential conflicting definitions.

template<class Policy, class InIt, class T, class BinOp>
inline std::enable_if_t<std::is_execution_policy_v<Policy>, T>
reduce(Policy&&, InIt first, InIt last, T init, BinOp op)
{
    for (; first != last; ++first)
    {
        init = op(init, *first);
    }
    return init;
}

template<class Policy, class InIt, class T>
inline std::enable_if_t<std::is_execution_policy_v<Policy>, T>
reduce(Policy&&, InIt first, InIt last, T init)
{
    for (; first != last; ++first)
    {
        init = (init + *first);
    }
    return init;
}


template
<
    class Policy,
    class InIt1, class InIt2,
    class T,
    class BinOp1, class BinOp2
>
inline enable_if_t<is_execution_policy_v<Policy>, T>
transform_reduce
(
    Policy&&,
    InIt1 first1, InIt1 last1, InIt2 first2, T init, BinOp1 op1, BinOp2 op2
)
{
    for (; first1 != last1; ++first1, (void)++first2)
    {
        init = op1(init, op2(*first1, *first2));
    }
    return init;
}


template<class Policy, class InIt, class OutIt, class T>
inline enable_if_t<is_execution_policy_v<Policy>, OutIt>
exclusive_scan(Policy&&, InIt first, InIt last, OutIt out, T init)
{
    for (; first != last; ++first, (void)++out)
    {
        T tmp(init);
        init = init + *first;
        *out = tmp;
    }
    return out;
}

} // End namespace std

#endif  /* [Foam_execution_fallback] */

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Gets quite noisy...
// #ifdef Foam_execution_fallback
// #warning "using hand-rolled <execution> replacements"
// #endif

// Cleanup. It is currently not needed elsewhere.
#undef Foam_execution_fallback

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
