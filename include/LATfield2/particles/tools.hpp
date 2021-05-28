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

#endif

}

#endif
