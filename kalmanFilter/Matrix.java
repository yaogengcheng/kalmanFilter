public class Matrix {
    private final int mRows, mCols;
    private final double[] mData;  // mRows * mCols, row major

    public Matrix(int rows, int cols) {
        mRows = rows;
        mCols = cols;
        mData = new double[mRows * mCols];
    }

    public Matrix(Matrix other) {
        mRows = other.mRows;
        mCols = other.mCols;
        mData = new double[mRows * mCols];
        for (int i = 0; i < mData.length; ++i) mData[i] = other.mData[i];
    }

    public final void set(int row, int col, double value) {
        assert row < mRows;
        assert col < mCols;
        mData[row * mCols + col] = value;
    }

    public final void setAll(double... values) {
        assert values.length == mData.length;
        for (int i = 0; i < values.length; ++i) mData[i] = values[i];
    }

    public final double get(int row, int col) {
        assert row < mRows;
        assert col < mCols;
        return mData[row * mCols + col];
    }

    public void setIdentityMatrix() {
        assert mRows == mCols;
        for (int i = 0; i < mRows; ++i) {
            for (int j = 0; j < mCols; ++j) {
                if (i == j) {
                    set(i, j, 1.0);
                } else {
                    set(i, j, 0.0);
                }
            }
        }
    }


    public String toString() {
        StringBuilder b = new StringBuilder();
        for (int i = 0; i < mRows; ++i) {
            for (int j = 0; j < mCols; ++j) {
                if (j > 0) b.append(" ");
                b.append(String.format("%.2f", get(i, j)));
            }
            b.append("\n");
        }
        return b.toString();
    }

    /**
     * c = a + b
     *
     * @p a and @p b remain unchanged.
     */
    public static Matrix add(Matrix a, Matrix b) {
        assert a.mData.length == b.mData.length;
        Matrix c = new Matrix(a.mRows,a.mCols);
        for (int i = 0; i < a.mData.length; ++i) {
            c.mData[i] = a.mData[i] + b.mData[i];
        }
        return c;
    }

    /**
     * c = a - b
     *
     * @p a and @p b remain unchanged.
     */
    public static Matrix subtract(Matrix a, Matrix b) {
        assert a.mData.length == b.mData.length;
        Matrix c = new Matrix(a.mRows,a.mCols);
        for (int i = 0; i < a.mData.length; ++i) {
            c.mData[i] = a.mData[i] - b.mData[i];
        }
        return c;
    }

    /**
     * this = I - this
     */
    public void subtractFromIdentityMatrix() {
        assert mRows == mCols;
        for (int i = 0; i < mRows; ++i) {
            for (int j = 0; j < mCols; ++j) {
                if (i == j) {
                    set(i, j, 1.0 - get(i, j));
                } else {
                    set(i, j, 0.0 - get(i, j));
                }
            }
        }
    }

    /**
     * c = a * b
     *
     * @p a and @p b remain unchanged.
     */
    public static Matrix multiply(Matrix a, Matrix b) {
        assert(a.mCols == b.mRows);
        Matrix c = new Matrix(a.mRows,b.mCols);
        for (int i = 0; i < c.mRows; ++i) {
            for (int j = 0; j < c.mCols; ++j) {
                // Calculate element c.data[i][j] via a dot product of one row of a
                // with one column of b
                double value = 0.0;
                for (int k = 0; k < a.mCols; ++k) {
                    value += a.get(i, k) * b.get(k, j);
                }
                c.set(i, j, value);
            }
        }
        return c;
    }

    /**
     * c = a * transpose(b)
     *
     * @p a and @p b remain unchanged.
     */
    public final static Matrix multiplyByTransposeMatrix(Matrix a, Matrix b) {
        assert(a.mCols == b.mCols);
        Matrix c = new Matrix(a.mRows,b.mRows);
        for (int i = 0; i < c.mRows; ++i) {
            for (int j = 0; j < c.mCols; ++j) {
                // Calculate element c.data[i][j] via a dot product of one row of a
                // with one row of b
                double value = 0.0;
                for (int k = 0; k < a.mCols; ++k) {
                    value += a.get(i, k) * b.get(j, k);
                }
                c.set(i, j, value);
            }
        }
        return c;
    }

    /**
     * output = transpose(input)
     */
    public final static Matrix transpose(Matrix input) {
        Matrix output = new Matrix(input.mCols,input.mRows);
        for (int i = 0; i < input.mRows; ++i) {
            for (int j = 0; j < input.mCols; ++j) {
                output.set(j, i, input.get(i, j));
            }
        }
        return output;
    }

    public final boolean equals(Matrix b, double tolerance) {
        assert(mData.length == b.mData.length);
        for (int i = 0; i < mData.length; ++i) {
            if (Math.abs(mData[i] - b.mData[i]) > tolerance) return false;
        }
        return true;
    }

    /**
     * this = this * scalar
     */
    public final void scale(double scalar) {
        assert(scalar != 0.0);
        for (int i = 0; i < mData.length; ++i) {
            mData[i] *= scalar;
        }
    }

    public final void swapRows(int r1, int r2) {
        assert(r1 != r2);
        for (int col = 0; col < mCols; ++col) {
            final double tmp = get(r1, col);
            set(r1, col, get(r2, col));
            set(r2, col, tmp);
        }
    }

    public final void scaleRow(int row, double scalar) {
        assert(scalar != 0.0);
        for (int col = 0; col < mCols; ++col) {
            set(row, col, get(row, col) * scalar);
        }
    }

    /* Add scalar * row r2 to row r1. */
    public void shearRow(int row1, int row2, double scalar) {
        assert(row1 != row2);
        for (int col = 0; col < mCols; ++col) {
            set(row1, col, get(row1, col) + scalar * get(row2, col));
        }
    }

    /** Uses Gauss-Jordan elimination.
     The elimination procedure works by applying elementary row
     operations to our input matrix until the input matrix is reduced to
     the identity matrix.
     Simultaneously, we apply the same elementary row operations to a
     separate identity matrix to produce the inverse matrix.
     If this makes no sense, read wikipedia on Gauss-Jordan elimination.

     This is not the fastest way to invert matrices, so this is quite
     possibly the bottleneck.
     */
    public static Matrix invert(Matrix in) {
        assert(in.mRows == in.mCols);
        Matrix output = new Matrix(in.mRows, in.mCols);
        output.setIdentityMatrix();
        Matrix input = new Matrix(in);

		  /* Convert input to the identity matrix via elementary row operations.
		     The ith pass through this loop turns the element at i,i to a 1
		     and turns all other elements in column i to a 0. */
        for (int i = 0; i < input.mRows; ++i) {
            if (input.get(i, i) == 0.0) {
                /* We must swap mRows to get a nonzero diagonal element. */
                int r;
                for (r = i + 1; r < input.mRows; ++r) {
                    if (input.get(r, i) != 0.0) {
                        break;
                    }
                }
                if (r == input.mRows) {
                    // Every remaining element in this column is zero, so this
                    // matrix cannot be inverted.
                    return null;
                }
                input.swapRows(i, r);
                output.swapRows(i, r);
            }

			  /* Scale this row to ensure a 1 along the diagonal.
		       We might need to worry about overflow from a huge scalar here. */
            double scalar = 1.0 / input.get(i, i);
            input.scaleRow(i, scalar);
            output.scaleRow(i, scalar);

            /* Zero out the other elements in this column. */
            for (int j = 0; j < input.mRows; ++j) {
                if (i == j) {
                    continue;
                }
                double shear_needed = -input.get(j, i);
                input.shearRow(j, i, shear_needed);
                output.shearRow(j, i, shear_needed);
            }
        }
        return output;
    }

    public void Zero() {
        for (double mDatum : mData) {
            mDatum = 0;
        }
    }
}