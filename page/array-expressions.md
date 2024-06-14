## Array expressions

Symbolics can represent variables that represent multi-dimensional arrays and operations on them. The operations are represented compactly and in constant space independent of the size of the array they result in. To create a variable that rerpesents an array, we still use the `@variables` macro but encode the size:

\repl{
@variables A[1:4, 1:5] b[1:5] c[1:4]
}

Operations on these variables create array expressions:

\repl{
A * b
A \ c
}

Note that this is different from an array of symbolic expressions. Native Julia `Array`s will work with most generic functions written intended for arrays (arrays of Real specifically if the expressions are all numeric).
