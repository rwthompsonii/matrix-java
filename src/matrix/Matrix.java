/*
 * Copyright (C) 2015 Ron
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
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package matrix;

import org.apache.commons.math3.util.Precision;

/**
 *
 * @author ron
 */
public class Matrix implements MatrixConstants {

    protected int rows;
    protected int columns;
    protected double[][] matrix;

    public int getRows() {
        return rows;
    }

    public int getColumns() {
        return columns;
    }

    public double[][] getMatrix() {
        return matrix;
    }

    public Matrix() {
        this.matrix = IdentityMatrix(DEFAULT_SIZE).getMatrix();
        this.rows = DEFAULT_SIZE;
        this.columns = DEFAULT_SIZE;
    }

    public Matrix(Matrix matrix) {
        this.matrix = matrix.getMatrix();
        this.rows = matrix.getRows();
        this.columns = matrix.getColumns();
    }

    public Matrix(double[][] matrix) {

        this.matrix = matrix;
        this.rows = matrix.length;
        this.columns = matrix[0].length;
    }

    public Matrix(double[] vector, boolean rowOrColumnVector) {
        if (rowOrColumnVector) {
            this.rows = 1;
            this.columns = vector.length;
            this.matrix = new double[rows][columns];
            for (int j = 0; j < columns; ++j) {
                matrix[0][j] = vector[j];
            }
        } else {
            this.columns = 1;
            this.rows = vector.length;
            this.matrix = new double[rows][columns];
            for(int i = 0; i < rows; ++i) {
                matrix[i][0] = vector[i];
            }
        }
    }

    /**
     * a convenience constructor provided for integer[][] inputs
     *
     * @param matrix
     */
    public Matrix(int[][] matrix) {

        this.rows = matrix.length;
        this.columns = matrix[0].length;

        double[][] convertedMatrix = new double[rows][columns];

        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < columns; ++j) {
                convertedMatrix[i][j] = (double) matrix[i][j];
            }
        }
        this.matrix = convertedMatrix;
    }

    /**
     * A convenience constructor provided for float[][] inputs
     *
     * @param matrix
     */
    public Matrix(float[][] matrix) {
        this.rows = matrix.length;
        this.columns = matrix[0].length;

        double[][] convertedMatrix = new double[rows][columns];

        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < columns; ++j) {
                convertedMatrix[i][j] = (double) matrix[i][j];
            }
        }
        this.matrix = convertedMatrix;
    }

    public static Matrix IdentityMatrix(int size) {
        assert (size > 1);

        double[][] identity = new double[size][size];

        //cannot use Arrays.fill on multidimensional arrays.... 
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                if (i == j) {
                    identity[i][j] = 1;
                } else {
                    identity[i][j] = 0;
                }
            }
        }

        Matrix Identity = new Matrix(identity);

        return Identity;
    }

    //note that equals does NOT compare with equality operator, it uses MatrixConstants.EPSILON for "close enough"
    public boolean equals(Matrix other) {
        if (other.getClass() != getClass() || other.getColumns() != getColumns() || other.getRows() != getRows()) {
            return false;
        } else {
            double[][] otherMatrix = other.getMatrix();
            for (int i = 0; i < rows; ++i) {
                for (int j = 0; j < columns; ++j) {
                    if (Math.abs(otherMatrix[i][j] - matrix[i][j]) > EPSILON) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    //this method is useful if you want to see the nasty decimals the default
    //toString() method hides
    public String rawToString() {
        StringBuilder result = new StringBuilder(rows * columns * 3);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < columns; ++j) {
                result.append(matrix[i][j]);
                if (j < columns - 1) {
                    result.append(",\t");
                }
            }
            result.append(";\n");
        }
        return result.toString();
    }

    @Override
    public String toString() {
        StringBuilder result = new StringBuilder(rows * columns * 3);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < columns; ++j) {
                //ugly floating point hack incoming to prevent 1.00000000004 and 8E-16 
                // and the like from crapping all over my pretty matrices
                result.append(Math.floor((matrix[i][j]) * 1000000000 + 0.5) / 1000000000);
                if (j < columns - 1) {
                    result.append(",\t");
                }
            }
            result.append(";\n");
        }
        return result.toString();
    }

    public String toCopyableString() {
        StringBuilder result = new StringBuilder(rows * columns * 3);
        result.append("{");
        for (int i = 0; i < rows; ++i) {
            result.append("{");
            for (int j = 0; j < columns; ++j) {
                result.append(Precision.round(matrix[i][j], 3));
                if (j < columns - 1) {
                    result.append(",");
                }
            }
            result.append("}");
            if (i < rows - 1) {
                result.append(",");
            }
        }
        result.append("}");
        return result.toString();
    }

    /**
     * function mult() is a convenience wrapper method used when you don't want
     * to create Matrix objects yourself to multiply two multidimensional
     * arrays. the objects are internally created to carry out the
     * multiplication, though, so it's just a wrapper in the true sense of the
     * word
     *
     * @param a - the left side operand of the matrix multiplication operation -
     * note matrix multiplication does not normally commute => in general => a*b
     * != b*a -- you were warned
     * @param b - the right side operand of the matrix multiplication operation
     * @return - a double[a.rows][b.columns] result array of the results of the
     * multiplication
     * @throws DimensionMismatchException - note that a.columns != b.rows will
     * result in this exception being thrown
     */
    public static double[][] mult(double[][] a, double[][] b) throws DimensionMismatchException {
        Matrix A = new Matrix(a);
        Matrix B = new Matrix(b);

        return Matrix.mult(A, B).getMatrix();
    }

    public Matrix mult(Matrix B) throws DimensionMismatchException {
        return mult(this, B);
    }

    public static Matrix mult(Matrix A, Matrix B) throws DimensionMismatchException {
        //you can think of this operation as A * B where A,B are matrices
        //remember that Matrix multiplication normally does NOT commute

        if (A.getColumns() != B.getRows()) {
            throw new DimensionMismatchException("The two matrices have invalid dimensions for multiplication.");
        }

        double[][] result = new double[A.getRows()][B.getColumns()];

        int aRows = A.getRows();
        int aCols = A.getColumns();
        int bCols = B.getColumns();

        double[][] a = A.getMatrix();
        double[][] b = B.getMatrix();

        //iterate over the resulting array filling in the values 
        //from the two operands according to matrix mult rules
        for (int i = 0; i < aRows; ++i) {
            for (int j = 0; j < bCols; ++j) {
                for (int k = 0; k < aCols/*could also use bRows here*/; ++k) {
                    result[i][j] += a[i][k] * b[k][j];
                }
            }
        }

        return new Matrix(result);
    }

    public static Matrix transpose(Matrix A) {
        double[][] result = new double[A.getColumns()][A.getRows()];

        double[][] a = A.getMatrix();

        for (int i = 0; i < A.getRows(); ++i) {
            for (int j = 0; j < A.getColumns(); ++j) {
                result[i][j] = a[j][i];
            }
        }

        return new Matrix(result);
    }

    public Matrix transpose() {
        return transpose(this);
    }

    public static Matrix scale(Matrix A, double v) {

        double[][] a = A.getMatrix().clone(); //using clone to avoid modifying A

        for (int i = 0; i < A.getRows(); ++i) {
            for (int j = 0; j < A.getColumns(); ++j) {
                a[i][j] *= v;
            }
        }

        return new Matrix(a);
    }

    public Matrix scale(double v) {
        return scale(this, v);
    }

    public static Matrix add(Matrix A, Matrix B) throws DimensionMismatchException {
        if (A.getRows() != B.getRows() || A.getColumns() != B.getColumns()) {
            throw new DimensionMismatchException("This operation can only be performed on matrices of equal dimensions.");
        }

        double[][] result = new double[A.getRows()][A.getColumns()];

        for (int i = 0; i < A.getRows(); ++i) {
            for (int j = 0; j < A.getColumns(); ++j) {
                result[i][j] = A.getMatrix()[i][j] + B.getMatrix()[i][j];
            }
        }

        return new Matrix(result);
    }

    public Matrix add(Matrix B) throws DimensionMismatchException {
        return add(this, B);
    }

    public static Matrix subtract(Matrix A, Matrix B) throws DimensionMismatchException {
        Matrix C = scale(B, -1);
        return add(A, C);
    }

    public Matrix subtract(Matrix B) throws DimensionMismatchException {
        return subtract(this, B);
    }

    public static Matrix zeroes(int rows, int columns) {
        double[][] zeroMatrix = new double[rows][columns];

        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < columns; ++j) {
                zeroMatrix[i][j] = 0;
            }
        }
        return new Matrix(zeroMatrix);
    }

    public static Matrix random(int rows, int columns, double max, double min) {
        double[][] randomMatrix = new double[rows][columns];

        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < columns; ++j) {
                randomMatrix[i][j] = (Math.random() * (max - min) + min);
            }
        }

        return new Matrix(randomMatrix);
    }
}
