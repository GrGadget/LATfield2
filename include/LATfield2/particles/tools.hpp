#ifndef PARTICLES_TOOLS_HPP
#define PARTICLES_TOOLS_HPP

#include "LATfield2/macros.hpp"
#include "LATfield2/Imag.hpp"
#include <iostream>
#include <type_traits>
#if __cplusplus < 201700L
namespace std
{
    template<typename ...> using void_t = void;
};
#endif

namespace LATfield2
{

#ifndef DOXYGEN_SHOULD_SKIP_THIS
LATfield2::Real get_lattice_resolution(const int npts[3],const Real boxSize[3]);

#define CREATE_MEMBER_DETECTOR(X)                \
template<typename T, typename = std::void_t<> >  \
struct has_##X : std::false_type                 \
{};                                              \
template<typename T>                             \
struct has_##X<T, std::void_t<decltype(T::X)> >  \
    : std::true_type                             \
{};                                              \

#define CREATE_MEMBER_DETECTOR_MAXI(X)           \
template<class P, class I>                       \
typename std::enable_if< has_##X<P>::value       \
    , double>::type                              \
get_##X( const P& p,const I&)                    \
{ return p.X; }                                  \
template<class P, class I>                       \
typename std::enable_if<                         \
    has_##X<P>::value==false                     \
    and has_##X<I>::value ,double>::type         \
get_##X(const P&,const I& i)                     \
{  return i.X; }                                 \
template<class P, class I>                       \
typename std::enable_if< has_##X<P>::value       \
    ,double&>::type                              \
get_##X(P& p,const I&)                           \
{ return p.X;}                                   \
template<class P, class I>                       \
typename std::enable_if<                         \
    has_##X<P>::value==false                     \
    and has_##X<I>::value ,double&>::type        \
get_##X(const P&,I& i)                           \
{ return i.X; }                                  \
template<class I>                                \
typename std::enable_if<                         \
    has_##X<I>::value==true ,double>::type       \
get_##X(const I& i)                              \
{  return i.X; }                                 \
template<class I>                                \
typename std::enable_if<                         \
    has_##X<I>::value==false ,double>::type      \
get_##X(const I& i)                              \
{  return -1.0; }                                \

//#define CREATE_MEMBER_DETECTOR(X)            \
//  template<typename T> struct has_##X {      \
//    struct Fallback {int X; };               \
//    struct Derived : T, Fallback { };        \
//    template<typename C, C> struct ChT;      \
//    template<typename C> static char (&f(ChT<int Fallback::*, &C::X>*))[1]; \
//    template<typename C> static char (&f(...))[2]; \
//    static bool const value = sizeof(f<Derived>(0)) == 2; \
//  };


// #define CREATE_MEMBER_DETECTOR_MAXI(X)          \
//     template<typename T> struct has_maxi_##X {      \
//     struct Fallback {int X; };               \
//     struct Derived : T, Fallback { };        \
//     template<typename C, C> struct ChT;      \
//     template<typename C,typename CC> int (f(ChT<int Fallback::*, &C::X>*)){return -1;}; \
//     template<typename C,typename CC> int (f(...)){return offsetof(T,X);}; \
//     int gos(){return f<Derived,T>(0);} \
// };

#endif

}

#endif
