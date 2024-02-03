#include "numc.h"
#include <structmember.h>

PyTypeObject Matrix61cType;

/* Helper functions for initalization of matrices and vectors */

/*
 * Return a tuple given rows and cols
 */
PyObject *get_shape(int rows, int cols) {
  if (rows == 1 || cols == 1) {
    return PyTuple_Pack(1, PyLong_FromLong(rows * cols));
  } else {
    return PyTuple_Pack(2, PyLong_FromLong(rows), PyLong_FromLong(cols));
  }
}
/*
 * Matrix(rows, cols, low, high). Fill a matrix random double values
 */
int init_rand(PyObject *self, int rows, int cols, unsigned int seed, double low,
              double high) {
    matrix *new_mat;
    int alloc_failed = allocate_matrix(&new_mat, rows, cols);
    if (alloc_failed) return alloc_failed;
    rand_matrix(new_mat, seed, low, high);
    ((Matrix61c *)self)->mat = new_mat;
    ((Matrix61c *)self)->shape = get_shape(new_mat->rows, new_mat->cols);
    return 0;
}

/*
 * Matrix(rows, cols, val). Fill a matrix of dimension rows * cols with val
 */
int init_fill(PyObject *self, int rows, int cols, double val) {
    matrix *new_mat;
    int alloc_failed = allocate_matrix(&new_mat, rows, cols);
    if (alloc_failed)
        return alloc_failed;
    else {
        fill_matrix(new_mat, val);
        ((Matrix61c *)self)->mat = new_mat;
        ((Matrix61c *)self)->shape = get_shape(new_mat->rows, new_mat->cols);
    }
    return 0;
}

/*
 * Matrix(rows, cols, 1d_list). Fill a matrix with dimension rows * cols with 1d_list values
 */
int init_1d(PyObject *self, int rows, int cols, PyObject *lst) {
    if (rows * cols != PyList_Size(lst)) {
        PyErr_SetString(PyExc_ValueError, "Incorrect number of elements in list");
        return -1;
    }
    matrix *new_mat;
    int alloc_failed = allocate_matrix(&new_mat, rows, cols);
    if (alloc_failed) return alloc_failed;
    int count = 0;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            set(new_mat, i, j, PyFloat_AsDouble(PyList_GetItem(lst, count)));
            count++;
        }
    }
    ((Matrix61c *)self)->mat = new_mat;
    ((Matrix61c *)self)->shape = get_shape(new_mat->rows, new_mat->cols);
    return 0;
}

/*
 * Matrix(2d_list). Fill a matrix with dimension len(2d_list) * len(2d_list[0])
 */
int init_2d(PyObject *self, PyObject *lst) {
    int rows = PyList_Size(lst);
    if (rows == 0) {
        PyErr_SetString(PyExc_ValueError,
                        "Cannot initialize numc.Matrix with an empty list");
        return -1;
    }
    int cols;
    if (!PyList_Check(PyList_GetItem(lst, 0))) {
        PyErr_SetString(PyExc_ValueError, "List values not valid");
        return -1;
    } else {
        cols = PyList_Size(PyList_GetItem(lst, 0));
    }
    for (int i = 0; i < rows; i++) {
        if (!PyList_Check(PyList_GetItem(lst, i)) ||
                PyList_Size(PyList_GetItem(lst, i)) != cols) {
            PyErr_SetString(PyExc_ValueError, "List values not valid");
            return -1;
        }
    }
    matrix *new_mat;
    int alloc_failed = allocate_matrix(&new_mat, rows, cols);
    if (alloc_failed) return alloc_failed;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            set(new_mat, i, j,
                PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(lst, i), j)));
        }
    }
    ((Matrix61c *)self)->mat = new_mat;
    ((Matrix61c *)self)->shape = get_shape(new_mat->rows, new_mat->cols);
    return 0;
}

/*
 * This deallocation function is called when reference count is 0
 */
void Matrix61c_dealloc(Matrix61c *self) {
    deallocate_matrix(self->mat);
    Py_TYPE(self)->tp_free(self);
}

/* For immutable types all initializations should take place in tp_new */
PyObject *Matrix61c_new(PyTypeObject *type, PyObject *args,
                        PyObject *kwds) {
    /* size of allocated memory is tp_basicsize + nitems*tp_itemsize*/
    Matrix61c *self = (Matrix61c *)type->tp_alloc(type, 0);
    return (PyObject *)self;
}

/*
 * This matrix61c type is mutable, so needs init function. Return 0 on success otherwise -1
 */
int Matrix61c_init(PyObject *self, PyObject *args, PyObject *kwds) {
    /* Generate random matrices */
    if (kwds != NULL) {
        PyObject *rand = PyDict_GetItemString(kwds, "rand");
        if (!rand) {
            PyErr_SetString(PyExc_TypeError, "Invalid arguments");
            return -1;
        }
        if (!PyBool_Check(rand)) {
            PyErr_SetString(PyExc_TypeError, "Invalid arguments");
            return -1;
        }
        if (rand != Py_True) {
            PyErr_SetString(PyExc_TypeError, "Invalid arguments");
            return -1;
        }

        PyObject *low = PyDict_GetItemString(kwds, "low");
        PyObject *high = PyDict_GetItemString(kwds, "high");
        PyObject *seed = PyDict_GetItemString(kwds, "seed");
        double double_low = 0;
        double double_high = 1;
        unsigned int unsigned_seed = 0;

        if (low) {
            if (PyFloat_Check(low)) {
                double_low = PyFloat_AsDouble(low);
            } else if (PyLong_Check(low)) {
                double_low = PyLong_AsLong(low);
            }
        }

        if (high) {
            if (PyFloat_Check(high)) {
                double_high = PyFloat_AsDouble(high);
            } else if (PyLong_Check(high)) {
                double_high = PyLong_AsLong(high);
            }
        }

        if (double_low >= double_high) {
            PyErr_SetString(PyExc_TypeError, "Invalid arguments");
            return -1;
        }

        // Set seed if argument exists
        if (seed) {
            if (PyLong_Check(seed)) {
                unsigned_seed = PyLong_AsUnsignedLong(seed);
            }
        }

        PyObject *rows = NULL;
        PyObject *cols = NULL;
        if (PyArg_UnpackTuple(args, "args", 2, 2, &rows, &cols)) {
            if (rows && cols && PyLong_Check(rows) && PyLong_Check(cols)) {
                return init_rand(self, PyLong_AsLong(rows), PyLong_AsLong(cols), unsigned_seed, double_low,
                                 double_high);
            }
        } else {
            PyErr_SetString(PyExc_TypeError, "Invalid arguments");
            return -1;
        }
    }
    PyObject *arg1 = NULL;
    PyObject *arg2 = NULL;
    PyObject *arg3 = NULL;
    if (PyArg_UnpackTuple(args, "args", 1, 3, &arg1, &arg2, &arg3)) {
        /* arguments are (rows, cols, val) */
        if (arg1 && arg2 && arg3 && PyLong_Check(arg1) && PyLong_Check(arg2) && (PyLong_Check(arg3)
                || PyFloat_Check(arg3))) {
            if (PyLong_Check(arg3)) {
                return init_fill(self, PyLong_AsLong(arg1), PyLong_AsLong(arg2), PyLong_AsLong(arg3));
            } else
                return init_fill(self, PyLong_AsLong(arg1), PyLong_AsLong(arg2), PyFloat_AsDouble(arg3));
        } else if (arg1 && arg2 && arg3 && PyLong_Check(arg1) && PyLong_Check(arg2) && PyList_Check(arg3)) {
            /* Matrix(rows, cols, 1D list) */
            return init_1d(self, PyLong_AsLong(arg1), PyLong_AsLong(arg2), arg3);
        } else if (arg1 && PyList_Check(arg1) && arg2 == NULL && arg3 == NULL) {
            /* Matrix(rows, cols, 1D list) */
            return init_2d(self, arg1);
        } else if (arg1 && arg2 && PyLong_Check(arg1) && PyLong_Check(arg2) && arg3 == NULL) {
            /* Matrix(rows, cols, 1D list) */
            return init_fill(self, PyLong_AsLong(arg1), PyLong_AsLong(arg2), 0);
        } else {
            PyErr_SetString(PyExc_TypeError, "Invalid arguments");
            return -1;
        }
    } else {
        PyErr_SetString(PyExc_TypeError, "Invalid arguments");
        return -1;
    }
}

/*
 * List of lists representations for matrices
 */
PyObject *Matrix61c_to_list(Matrix61c *self) {
    int rows = self->mat->rows;
    int cols = self->mat->cols;
    PyObject *py_lst = NULL;
    if (self->mat->is_1d) {  // If 1D matrix, print as a single list
        py_lst = PyList_New(rows * cols);
        int count = 0;
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                PyList_SetItem(py_lst, count, PyFloat_FromDouble(get(self->mat, i, j)));
                count++;
            }
        }
    } else {  // if 2D, print as nested list
        py_lst = PyList_New(rows);
        for (int i = 0; i < rows; i++) {
            PyList_SetItem(py_lst, i, PyList_New(cols));
            PyObject *curr_row = PyList_GetItem(py_lst, i);
            for (int j = 0; j < cols; j++) {
                PyList_SetItem(curr_row, j, PyFloat_FromDouble(get(self->mat, i, j)));
            }
        }
    }
    return py_lst;
}

PyObject *Matrix61c_class_to_list(Matrix61c *self, PyObject *args) {
    PyObject *mat = NULL;
    if (PyArg_UnpackTuple(args, "args", 1, 1, &mat)) {
        if (!PyObject_TypeCheck(mat, &Matrix61cType)) {
            PyErr_SetString(PyExc_TypeError, "Argument must of type numc.Matrix!");
            return NULL;
        }
        Matrix61c* mat61c = (Matrix61c*)mat;
        return Matrix61c_to_list(mat61c);
    } else {
        PyErr_SetString(PyExc_TypeError, "Invalid arguments");
        return NULL;
    }
}

/*
 * Add class methods
 */
PyMethodDef Matrix61c_class_methods[] = {
    {"to_list", (PyCFunction)Matrix61c_class_to_list, METH_VARARGS, "Returns a list representation of numc.Matrix"},
    {NULL, NULL, 0, NULL}
};

/*
 * Matrix61c string representation. For printing purposes.
 */
PyObject *Matrix61c_repr(PyObject *self) {
    PyObject *py_lst = Matrix61c_to_list((Matrix61c *)self);
    return PyObject_Repr(py_lst);
}

/* NUMBER METHODS */

/*
 * Add the second numc.Matrix (Matrix61c) object to the first one. The first operand is
 * self, and the second operand can be obtained by casting `args`.
 */
PyObject *Matrix61c_add(Matrix61c* self, PyObject* args) {
    /* TODO: YOUR CODE HERE */
    if (!PyObject_TypeCheck(args, &Matrix61cType)) {
        PyErr_SetString(PyExc_TypeError, "Argument must be of type numc.Matrix!");
        return Py_None;
    }
    if (self->mat->rows != ((Matrix61c*)args)->mat->rows || self->mat->cols != ((Matrix61c*)args)->mat->cols) {
        PyErr_SetString(PyExc_ValueError, "Matrices should have the same dimensions!");
        return Py_None;
    }

    // create a new Matrix object
    Matrix61c* result = (Matrix61c*) Matrix61c_new(&Matrix61cType, NULL, NULL);
    PyObject* mat_args = PyTuple_Pack(2, PyLong_FromLong(self->mat->rows), PyLong_FromLong(self->mat->cols));
    Matrix61c_init((PyObject*)result, mat_args, NULL);
    // operation add
    add_matrix(result->mat, self->mat, ((Matrix61c*)args)->mat);
    result->shape = get_shape(self->mat->rows, self->mat->cols);
    return result;
}

/*
 * Substract the second numc.Matrix (Matrix61c) object from the first one. The first operand is
 * self, and the second operand can be obtained by casting `args`.
 */
PyObject *Matrix61c_sub(Matrix61c* self, PyObject* args) {
    /* TODO: YOUR CODE HERE */
    if (!PyObject_TypeCheck(args, &Matrix61cType)) {
        PyErr_SetString(PyExc_TypeError, "Argument must be of type numc.Matrix!");
        return Py_None;
    }
    if (self->mat->rows != ((Matrix61c*)args)->mat->rows || self->mat->cols != ((Matrix61c*)args)->mat->cols) {
        PyErr_SetString(PyExc_ValueError, "Matrices should have the same dimensions!");
        return Py_None;
    }

    // create a new Matrix object
    Matrix61c* result = Matrix61c_new(&Matrix61cType, NULL, NULL);
    PyTupleObject* mat_args = PyTuple_Pack(2, PyLong_FromLong(self->mat->rows), PyLong_FromLong(self->mat->cols));
    Matrix61c_init(result, mat_args, NULL);
    // operation sub
    sub_matrix(result->mat, self->mat, ((Matrix61c*)args)->mat);
    result->shape = get_shape(self->mat->rows, self->mat->cols);
    return result;
}

/*
 * NOT element-wise multiplication. The first operand is self, and the second operand
 * can be obtained by casting `args`.
 */
PyObject *Matrix61c_multiply(Matrix61c* self, PyObject *args) {
    /* TODO: YOUR CODE HERE */
    if (!PyObject_TypeCheck(args, &Matrix61cType)) {
        PyErr_SetString(PyExc_TypeError, "Argument must be of type numc.Matrix!");
        return Py_None;
    }
    if (self->mat->cols != ((Matrix61c*)args)->mat->rows) {
        PyErr_SetString(PyExc_ValueError, "Number of former columns is not equal to number of latter rows!");
        return Py_None;
    }
    Matrix61c* result = Matrix61c_new(&Matrix61cType, NULL, NULL);
    PyTupleObject* mat_args = PyTuple_Pack(2, PyLong_FromLong(self->mat->rows), PyLong_FromLong(self->mat->cols));
    Matrix61c_init(result, mat_args, NULL);
    // operaton muliply
    mul_matrix(result->mat, self->mat, ((Matrix61c*)args)->mat);
    result->shape = get_shape(self->mat->rows, self->mat->cols);
    return result;
}

/*
 * Negates the given numc.Matrix.
 */
PyObject *Matrix61c_neg(Matrix61c* self) {
    /* TODO: YOUR CODE HERE */
    Matrix61c* result = Matrix61c_new(&Matrix61cType, NULL, NULL);
    PyTupleObject* mat_args = PyTuple_Pack(2, PyLong_FromLong(self->mat->rows), PyLong_FromLong(self->mat->cols));
    Matrix61c_init(result, mat_args, NULL);
    // operation neg
    neg_matrix(result->mat, self->mat);
    result->shape = get_shape(self->mat->rows, self->mat->cols);
    return result;
}

/*
 * Take the element-wise absolute value of this numc.Matrix.
 */
PyObject *Matrix61c_abs(Matrix61c *self) {
    /* TODO: YOUR CODE HERE */
    Matrix61c* result = Matrix61c_new(&Matrix61cType, NULL, NULL);
    PyTupleObject* mat_args = PyTuple_Pack(2, PyLong_FromLong(self->mat->rows), PyLong_FromLong(self->mat->cols));
    Matrix61c_init(result, mat_args, NULL);
    // operation abs
    abs_matrix(result->mat, self->mat);
    result->shape = get_shape(self->mat->rows, self->mat->cols);
    return result;
}

/*
 * Raise numc.Matrix (Matrix61c) to the `pow`th power. You can ignore the argument `optional`.
 */
PyObject *Matrix61c_pow(Matrix61c *self, PyObject *pow, PyObject *optional) {
    /* TODO: YOUR CODE HERE */
    if (!PyLong_Check(pow)) {
        PyErr_SetString(PyExc_TypeError, "Argument pow should be an integer!");
        return Py_None;
    }
    if (self->mat->rows != self->mat->cols) {
        PyErr_SetString(PyExc_ValueError, "Matrix should be square!");
        return Py_None;
    }
    if (PyLong_AsLong(pow) < 0) {
        PyErr_SetString(PyExc_ValueError, "Argument pow should be non-negative!");
        return Py_None;
    }
    Matrix61c* result = Matrix61c_new(&Matrix61cType, NULL, NULL);
    PyTupleObject* mat_args = PyTuple_Pack(2, PyLong_FromLong(self->mat->rows), PyLong_FromLong(self->mat->cols));
    Matrix61c_init(result, mat_args, NULL);
    // operation pow
    pow_matrix(result->mat, self->mat, PyLong_AsLong(pow));
    result->shape = get_shape(self->mat->rows, self->mat->cols);
    return result;
}

/*
 * Create a PyNumberMethods struct for overloading operators with all the number methods you have
 * define. You might find this link helpful: https://docs.python.org/3.6/c-api/typeobj.html
 */
PyNumberMethods Matrix61c_as_number = {
    /* TODO: YOUR CODE HERE */
    Matrix61c_add,
    Matrix61c_sub,
    Matrix61c_multiply,
    0,
    0,
    Matrix61c_pow,
    Matrix61c_neg,
    0,
    Matrix61c_abs,
};


/* INSTANCE METHODS */

/*
 * Given a numc.Matrix self, parse `args` to (int) row, (int) col, and (double/int) val.
 * Return None in Python (this is different from returning null).
 */
PyObject *Matrix61c_set_value(Matrix61c *self, PyObject* args) {
    /* TODO: YOUR CODE HERE */
    PyObject* row = NULL;
    PyObject* col = NULL;
    PyObject* val = NULL;
    if (!PyArg_UnpackTuple(args, "args", 3, 3, &row, &col, &val)) {
        PyErr_SetString(PyExc_TypeError, "The number of arguments parsed from args should be 3!");
        return Py_None;
    }
    if (!PyLong_Check(row) 
        || !PyLong_Check(col)) {
        PyErr_SetString(PyExc_TypeError, "Argument row and col should be integer!");
        return Py_None;
    }
    if (!PyLong_Check(val) 
        || !PyFloat_Check(val)) {
        PyErr_SetString(PyExc_TypeError, "The value should be float or integer!");
        return Py_None;
    }
    if (PyLong_AsLong(row) < 0 
        || PyLong_AsLong(row) >= self->mat->rows
        || PyLong_AsLong(col) < 0
        || PyLong_AsLong(col) >= self->mat->cols) {
        PyErr_SetString(PyExc_IndexError, "Accessing to Matrix[row, col] is out of range!");
        return Py_None;
    }
    // operation set
    if (PyLong_Check(val)) {
        set(self->mat, PyLong_AsLong(row), PyLong_AsLong(col), PyLong_AsLong(val));
    }
    else if (PyFloat_Check(val)) {
        set(self->mat, PyLong_AsLong(row), PyLong_AsLong(col), PyFloat_AsDouble(val));
    }
    else {
        PyErr_SetString(PyExc_TypeError, "The valu should be float or integer!");
    }
    return Py_None;
}

/*
 * Given a numc.Matrix `self`, parse `args` to (int) row and (int) col.
 * Return the value at the `row`th row and `col`th column, which is a Python
 * float/int.
 */
PyObject *Matrix61c_get_value(Matrix61c *self, PyObject* args) {
    /* TODO: YOUR CODE HERE */
    // PyObject* row = NULL;
    PyObject* row = Py_None;
    // PyObject* col = NULL;
    PyObject* col = Py_None;
    if (!PyArg_UnpackTuple(args, "args", 2, 2, &row, &col)) {
        PyErr_SetString(PyExc_TypeError, "The number of arguments parsed from args should be 2!");
        return Py_None;
    }
    if (!PyLong_Check(row) || !PyLong_Check(col)) {
        PyErr_SetString(PyExc_TypeError, "Argument row and col should be integer!");
        return Py_None;
    }
    if (PyLong_AsLong(row) < 0 
        || PyLong_AsLong(row) >= self->mat->rows
        || PyLong_AsLong(col) < 0
        || PyLong_AsLong(col) >= self->mat->cols) {
        PyErr_SetString(PyExc_IndexError, "Accessing to Matrix[row, col] is out of range!");
        return Py_None;
    }
    // operation get
    // PyObject* val = PyFloat_FromDouble(get(self->mat, PyLong_AsLong(row), PyLong_AsLong(col)));
    return PyFloat_FromDouble(get(self->mat, PyLong_AsLong(row), PyLong_AsLong(col)));
}

/*
 * Create an array of PyMethodDef structs to hold the instance methods.
 * Name the python function corresponding to Matrix61c_get_value as "get" and Matrix61c_set_value
 * as "set"
 * You might find this link helpful: https://docs.python.org/3.6/c-api/structures.html
 */
PyMethodDef Matrix61c_methods[] = {
    /* TODO: YOUR CODE HERE */
    {"set", (PyCFunction)Matrix61c_set_value, METH_VARARGS, "return None in Python"},
    {"get", (PyCFunction)Matrix61c_get_value, METH_VARARGS, "return the value at the 'row'th row and 'col'th column, which is a Python float/int"},
    {NULL, NULL, 0, NULL}
};

/* INDEXING */

/*
 * Given a numc.Matrix `self`, index into it with `key`. Return the indexed result.
 */

PyObject *Matrix61c_subscript(Matrix61c* self, PyObject* key) {
    /* TODO: YOUR CODE HERE */
    Matrix61c* subscript = (Matrix61c*) Matrix61c_new(&Matrix61cType, NULL, NULL);
    if (self->mat->is_1d) {
        // matrix is 1D
        // >>> b = nc.Matrix(1, 3) # b is a 1D matrix
        // or
        // >>> b = nc.Matrix(3, 1) # b is a 1D matrix
        if (PyLong_Check(key)) {
            // >>> b[0] # key is a single integer
            // 0.0

            int row = 1;
            int col = 1;
            if (self->mat->rows == 1) {
                col = PyLong_AsLong(key);
                if (col >= self->mat->cols) {
                    PyErr_SetString(PyExc_IndexError, "key is out of range!");
                    return subscript;
                }
            }
            else {
                row = PyLong_AsLong(key);
                if (row >= self->mat->rows) {
                    PyErr_SetString(PyExc_IndexError, "key is out of range!");
                    return subscript;
                }
            }
            allocate_matrix_ref(&subscript->mat, self->mat, 0, 0, row, col);
        }
        else if (PySlice_Check(key)) {
            // >>> b[0:2] # key is a single slice
            // [0.0, 0.0]

            Py_ssize_t start;
            Py_ssize_t stop;
            Py_ssize_t step;
            Py_ssize_t slicelength;

            int fail = 0;

            if (self->mat->rows == 1) {
                fail = PySlice_GetIndicesEx(key, self->mat->cols, &start, &stop, &step, &slicelength);
            }
            else {
                fail = PySlice_GetIndicesEx(key, self->mat->rows, &start, &stop, &step, &slicelength);
            }


            if (fail) {
                PyErr_SetString(PyExc_ValueError, "PySlice_GetINdicesEx() failed!");
                return subscript;
            }

            if (slicelength < 1) {
                PyErr_SetString(PyExc_ValueError, "the length of the slice < 1!");
                return subscript;
            }

            if (step != 1) {
                PyErr_SetString(PyExc_ValueError, "the step of the slice is not equal to 1!");
                return subscript;
            }

            if (self->mat->rows == 1) {
                allocate_matrix_ref(&subscript->mat, self->mat, 0, 0, 1, stop);
            }
            else {
                allocate_matrix_ref(&subscript->mat, self->mat, 0, 0, start, 1);
            }
        }
        else {
            // error handle
            // >>> b[0:1, 0:1] # This is invalid!
            // Traceback (most recent call last):
            // File "<stdin>", line 1, in <module>
            // TypeError: 1D matrices only support single slice!
            PyErr_SetString(PyExc_TypeError, "1D matrices only support single slice!");
            return subscript;
        }
    }
    else {
        // matrix is 2D
        // >>> a = nc.Matrix(3, 3)
        if (PyLong_Check(key)) {
            // >>> a[0] # key is a single number
            // [0.0, 0.0, 0.0]
            int index = PyLong_AsLong(key);
            if (index >= self->mat->rows) {
                PyErr_SetString(PyExc_IndexError, "key is out of range!");
                return subscript;
            }
            allocate_matrix_ref(&subscript->mat, self->mat, index, 0, 1, self->mat->cols);
        }
        else if (PySlice_Check(key)) {
            // todo: 
            // >>> a[0:2] # key is a single slice
            // [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
            Py_ssize_t start;
            Py_ssize_t stop;
            Py_ssize_t step;
            Py_ssize_t slicelength;
            int fail = PySlice_GetIndicesEx(key, self->mat->rows, &start, &stop, &step, &slicelength);

            if (fail) {
                PyErr_SetString(PyExc_ValueError, "key is not a useful slice!");
                return subscript;
            }

            if (start < 0 
                || start >= self->mat->rows 
                || stop < 1 
                || stop > self->mat->rows) 
            {
                PyErr_SetString(PyExc_IndexError, "index is out of range!");
                return subscript;
            }

            if (slicelength < 1) {
                PyErr_SetString(PyExc_ValueError, "the length of the slice < 1!");
                return subscript;
            }

            if (step != 1) {
                PyErr_SetString(PyExc_ValueError, "the step of the slice is not equal to 1!");
                return subscript;
            }

            allocate_matrix_ref(&subscript->mat, self->mat, start, 0, stop - start, self->mat->cols);

        }
        else if (PyTuple_Check(key)) {
            // todo: key is a tuple 
            if (PyTuple_GET_SIZE(key) != 2) {
                PyErr_SetString(PyExc_TypeError, "the length of tuple is not equal to 2!");
                return subscript;
            }

            PyObject* arg0 = PyTuple_GetItem(key, 0);
            PyObject* arg1 = PyTuple_GetItem(key, 1);


            if (PySlice_Check(arg0) && PySlice_Check(arg1)) {
                // >>> a[0:2, 0:2] # key is a tuple of two slices
                // [[0.0, 0.0], [0.0, 0.0]]
                Py_ssize_t row_start;
                Py_ssize_t row_end;
                Py_ssize_t row_slicelength;
                Py_ssize_t col_start;
                Py_ssize_t col_end;
                Py_ssize_t step;
                Py_ssize_t col_slicelength;

                int fail = PySlice_GetIndicesEx(arg0, self->mat->rows, &row_start, &row_end, &step, &row_slicelength);
                if (fail) {
                    PyErr_SetString(PyExc_ValueError, "PySlice_GetINdicesEx() failed!");
                    return subscript;
                }

                if (row_slicelength < 1) {
                    PyErr_SetString(PyExc_ValueError, "the length of the slice < 1!");
                    return subscript;
                }

                if (step != 1) {
                    PyErr_SetString(PyExc_ValueError, "the step of the slice is not equal to 1!");
                    return subscript;
                }

                fail = PySlice_GetIndicesEx(arg1, self->mat->cols, &col_start, &col_end, &step, &col_slicelength);
                if (fail) {
                    PyErr_SetString(PyExc_ValueError, "PySlice_GetINdicesEx() failed!");
                    return subscript;
                }

                if (col_slicelength < 1) {
                    PyErr_SetString(PyExc_ValueError, "the length of the slice < 1!");
                    return subscript;
                }

                if (step != 1) {
                    PyErr_SetString(PyExc_ValueError, "the step of the slice is not equal to 1!");
                    return subscript;
                }

                allocate_matrix_ref(&subscript->mat, self->mat, row_start, col_start, row_slicelength, col_slicelength);

            }
            else if (PySlice_Check(arg0) && PyLong_Check(arg1)) {
                // >>> a[0:2, 0] # key is a tuple of (slice, int)
                // [0.0, 0.0]
                Py_ssize_t row_start;
                Py_ssize_t row_end;
                Py_ssize_t row_slicelength;
                Py_ssize_t step;

                int fail = PySlice_GetIndicesEx(arg0, self->mat->rows, &row_start, &row_end, &step, &row_slicelength);
                if (fail) {
                    PyErr_SetString(PyExc_ValueError, "PySlice_GetINdicesEx() failed!");
                    return subscript;
                }

                if (row_slicelength < 1) {
                    PyErr_SetString(PyExc_ValueError, "the length of the slice < 1!");
                    return subscript;
                }

                if (step != 1) {
                    PyErr_SetString(PyExc_ValueError, "the step of the slice is not equal to 1!");
                    return subscript;
                }

                int col_start = PyLong_AsLong(arg1);

                if (col_start < 0 || col_start >= self->mat->cols) {
                    PyErr_SetString(PyExc_IndexError, "index is out of range!");
                    return subscript;
                }

                allocate_matrix_ref(&subscript->mat, self->mat, row_start, col_start, row_slicelength, 1);


            }
            else if (PyLong_Check(arg0) && PySlice_Check(arg1)) {
                // >>> a[0, 0:2] # key is a tuple of (int, slice)
                // [0.0, 0.0]
                int row_start = PyLong_AsLong(arg0);

                if (row_start < 0 || row_start >= self->mat->rows) {
                    PyErr_SetString(PyExc_IndexError, "index is out of range!");
                    return subscript;
                }

                Py_ssize_t col_start;
                Py_ssize_t col_end;
                Py_ssize_t step;
                Py_ssize_t col_slicelength;

                int fail = PySlice_GetIndicesEx(arg1, self->mat->cols, &col_start, &col_end, &step, &col_slicelength);
                if (fail) {
                    PyErr_SetString(PyExc_ValueError, "PySlice_GetINdicesEx() failed!");
                    return subscript;
                }

                if (col_slicelength < 1) {
                    PyErr_SetString(PyExc_ValueError, "the length of the slice < 1!");
                    return subscript;
                }

                if (step != 1) {
                    PyErr_SetString(PyExc_ValueError, "the step of the slice is not equal to 1!");
                    return subscript;
                }

                allocate_matrix_ref(&subscript->mat, self->mat, row_start, col_start, 1, col_slicelength);
            }
            else if (PyLong_Check(arg0) && PyLong_Check(arg1)) {
                // >>> a[0, 0] # key is a tuple of (int, int)
                // 0.0
                int row_start = PyLong_AsLong(arg0);

                if (row_start < 0 || row_start >= self->mat->rows) {
                    PyErr_SetString(PyExc_IndexError, "index is out of range!");
                    return subscript;
                }

                int col_start = PyLong_AsLong(arg1);

                if (col_start < 0 || col_start >= self->mat->cols) {
                    PyErr_SetString(PyExc_IndexError, "index is out of range!");
                    return subscript;
                }

                allocate_matrix_ref(&subscript->mat, self->mat, row_start, col_start, 1, 1);
            }
            else {
                // error handle
                PyErr_SetString(PyExc_TypeError, "the tuple is not of slices/ints!");
                return subscript;
            }
        }
        else {
            PyErr_SetString(PyExc_TypeError, "key is not an integer, a slice, or a length-2 tuple of slices/ints.!");
            return subscript;
        }
    }
    subscript->shape = get_shape(subscript->mat->rows, subscript->mat->cols);
    return subscript;
}

/*
 * Given a numc.Matrix `self`, index into it with `key`, and set the indexed result to `v`.
 */
int Matrix61c_set_subscript(Matrix61c* self, PyObject *key, PyObject *v) {
    /* TODO: YOUR CODE HERE */
    // get subscript matrix of the original
    Matrix61c* submat = (Matrix61c*) Matrix61c_subscript(self, key);
    if (submat->shape == Py_None) {
        return -1;
    }

    if (PyTuple_GET_SIZE(submat->shape) == 1) {
        // 1D submatrix
        if (PyLong_AsLong(PyTuple_GetItem(submat->shape, 0)) == 1) {
            // v should be a float or int
            if (PyLong_Check(v)) {
                submat->mat->data[0][0] = PyLong_AsLong(v);
            }
            else if (PyFloat_Check(v)) {
                submat->mat->data[0][0] = PyFloat_AsDouble(v);
            }
            else {
                PyErr_SetString(PyExc_TypeError, "v is not a float or int!");
                return -1;
            }
        }
        else {
            // v should be a list
            if (PyList_Check(v)) {
                PyErr_SetString(PyExc_TypeError, "v is not a list!");
                return -1;
            }

            if (PyList_GET_SIZE(v) != submat->mat->cols) {
                PyErr_SetString(PyExc_ValueError, "list has wrong size!");
                return -1;
            }

            for (int col = 0; col < submat->mat->cols; ++col) {
                PyObject* val = PyList_GetItem(v, col);
                if (PyLong_Check(val)) {
                    submat->mat->data[0][col] = PyLong_AsLong(val);
                }
                else if (PyFloat_Check(val)) {
                    submat->mat->data[0][col] = PyFloat_AsDouble(val);
                }
                else {
                    PyErr_SetString(PyExc_TypeError, "list has wrong size!");
                    return -1;
                }
            }
        }
    }
    else {
        // 2D submatrix
        if (!PyList_Check(v)) {
            PyErr_SetString(PyExc_TypeError, "v is not a list!");
            return -1;
        }

        if (PyList_GET_SIZE(v) != submat->mat->rows) {
            PyErr_SetString(PyExc_ValueError, "lis has wrong size!");
            return -1;
        }

        for (int row = 0; row < submat->mat->rows; ++row) {
            PyObject* lst = PyList_GetItem(v, row);
            if (!PyList_Check(lst)) {
                PyErr_SetString(PyExc_TypeError, "list is required!");
                return -1;
            }
            if (PyList_GET_SIZE(lst) != submat->mat->cols) {
                PyErr_SetString(PyExc_ValueError, "list has wrong size!");
                return -1;
            }
            for (int col = 0; col < submat->mat->cols; ++col) {
                PyObject* val = PyList_GetItem(lst, col);
                if (PyLong_Check(val)) {
                    submat->mat->data[row][col] = PyLong_AsLong(val);
                }
                else if (PyFloat_Check(val)) {
                    submat->mat->data[row][col] = PyFloat_AsDouble(val);
                }
                else {
                    PyErr_SetString(PyExc_TypeError, "val is not a float or int!");
                    return -1;
                }
            }
        }
    }
    return 0;
}

PyMappingMethods Matrix61c_mapping = {
    NULL,
    (binaryfunc) Matrix61c_subscript,
    (objobjargproc) Matrix61c_set_subscript,
};

/* INSTANCE ATTRIBUTES*/
PyMemberDef Matrix61c_members[] = {
    {
        "shape", T_OBJECT_EX, offsetof(Matrix61c, shape), 0,
        "(rows, cols)"
    },
    {NULL}  /* Sentinel */
};

PyTypeObject Matrix61cType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "numc.Matrix",
    .tp_basicsize = sizeof(Matrix61c),
    .tp_dealloc = (destructor)Matrix61c_dealloc,
    .tp_repr = (reprfunc)Matrix61c_repr,
    .tp_as_number = &Matrix61c_as_number,
    .tp_flags = Py_TPFLAGS_DEFAULT |
    Py_TPFLAGS_BASETYPE,
    .tp_doc = "numc.Matrix objects",
    .tp_methods = Matrix61c_methods,
    .tp_members = Matrix61c_members,
    .tp_as_mapping = &Matrix61c_mapping,
    .tp_init = (initproc)Matrix61c_init,
    .tp_new = Matrix61c_new
};


struct PyModuleDef numcmodule = {
    PyModuleDef_HEAD_INIT,
    "numc",
    "Numc matrix operations",
    -1,
    Matrix61c_class_methods
};

/* Initialize the numc module */
PyMODINIT_FUNC PyInit_numc(void) {
    PyObject* m;

    if (PyType_Ready(&Matrix61cType) < 0)
        return NULL;

    m = PyModule_Create(&numcmodule);
    if (m == NULL)
        return NULL;

    Py_INCREF(&Matrix61cType);
    PyModule_AddObject(m, "Matrix", (PyObject *)&Matrix61cType);
    printf("CS61C Fall 2020 Project 4: numc imported!\n");
    fflush(stdout);
    return m;
}