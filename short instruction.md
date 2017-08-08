
# Very short instruction

Objects are initialized directly (returns integers) using the init* functions, or indirectly using createObject or make.bspline, which return objects with class functionObject.

Evaluation is through the 'eval' functions: eval.fct and evl.deriv take functionObjects as arguments, while eval.direct and eval.deriv.direct take integers (memory addresses) as arguments,
check_validity (default = TRUE) checks if address points to a valid object before execution.




