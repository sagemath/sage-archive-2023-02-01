


#if !defined(MAX_MIN_H)
#define MAX_MIN_H


// Doubly templated version of max and min. =)

template<class Type1, class Type2>
inline
Type1 max(Type1 x, Type2 y) {
  if (x > y)
    return x;
  else
    return (Type1) y;
}


template<class Type1, class Type2>
inline
Type1 min(Type1 x, Type2 y) {
  if (x < y)
    return x;
  else
    return (Type1) y;
}


#endif
