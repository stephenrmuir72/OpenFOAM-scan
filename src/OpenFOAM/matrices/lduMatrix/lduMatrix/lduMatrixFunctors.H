#ifndef lduMatrixFunctors_H
#define lduMatrixFunctors_H

namespace Foam
{

#define GENERATE_UNARY_OPERATOR_FUNCTORS(Op,opFunc)                            \
template<class Type,class RType>                                               \
struct opFunc##UnaryOperatorFunctor                                            \
{                                                                              \
    RType operator()(const Type& t)                                            \
    {                                                                          \
        return Op t;                                                           \
    }                                                                          \
};                                                                             \

#define GENERATE_UNARY_FUNCTION_FUNCTORS(opFunc)                               \
template<class Type,class RType>                                               \
struct opFunc##UnaryFunctionFunctor                                            \
{                                                                              \
    RType operator()(const Type& t) const                                      \
    {                                                                          \
        return opFunc(t);                                                      \
    }                                                                          \
};                                                                             \

GENERATE_UNARY_OPERATOR_FUNCTORS(-, negate)
GENERATE_UNARY_OPERATOR_FUNCTORS(+, unity)
GENERATE_UNARY_FUNCTION_FUNCTORS(mag)

}
#endif
