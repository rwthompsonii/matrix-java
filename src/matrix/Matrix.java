/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package matrix;

/**
 *
 * @author ron
 */
public class Matrix {

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
        this.matrix = IdentityMatrix(2).getMatrix();
        this.rows = 2;
        this.columns = 2;
    }

    public Matrix(Matrix matrix) {
        this.matrix = matrix.getMatrix();
        this.rows = matrix.getRows();
        this.columns = matrix.getColumns();
    }

    public Matrix(double[][] matrix, int rows, int columns) {
        assert (rows == matrix.length);
        assert (columns == matrix[0].length);

        this.matrix = matrix;
        this.rows = rows;
        this.columns = columns;
    }

    /**
     * a convenience constructor provided for integer[][] inputs
     *
     * @param matrix
     * @param rows
     * @param columns
     */
    public Matrix(int[][] matrix, int rows, int columns) {
        assert (rows == matrix.length);
        assert (columns == matrix[0].length);

        double[][] convertedMatrix = new double[rows][columns];
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < columns; ++j) {
                convertedMatrix[i][j] = (double) matrix[i][j];
            }
        }

        this.matrix = convertedMatrix;
        this.rows = rows;
        this.columns = columns;
    }

    /**
     * A convenience constructor provided for float[][] inputs
     *
     * @param matrix
     * @param rows
     * @param columns
     */
    public Matrix(float[][] matrix, int rows, int columns) {
        assert (rows == matrix.length);
        assert (columns == matrix[0].length);

        double[][] convertedMatrix = new double[rows][columns];
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < columns; ++j) {
                convertedMatrix[i][j] = (double) matrix[i][j];
            }
        }

        this.matrix = convertedMatrix;
        this.rows = rows;
        this.columns = columns;
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

        Matrix Identity = new Matrix(identity, size, size);

        return Identity;
    }

    public boolean equals(Matrix other) {
        if (other.getClass() != getClass() || other.getColumns() != getColumns() || other.getRows() != getRows()) {
            return false;
        } else {
            double[][] otherMatrix = other.getMatrix();
            for (int i = 0; i < rows; ++i) {
                for (int j = 0; j < columns; ++j) {
                    if (otherMatrix[i][j] != matrix[i][j]) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    @Override
    public String toString() {
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

        //iterate over the resulting array filling in the values from the two operands according to matrix mult rules
        for (int i = 0; i < aRows; ++i) {
            for (int j = 0; j < bCols; ++j) {
                for (int k = 0; k < aCols/*could also use bRows here*/; ++k) {
                    result[i][j] += a[i][k] * b[k][j];
                }
                /*if (Math.abs(result[i][j]) < 0.00000001) {
                 result[i][j] = 0.0; //ugly hack to get around floating point operations having limited precision
                 }*/
                result[i][j] = Math.floor((result[i][j] * 1000000000 + 0.5) / 1000000000);
                //slightly uglier hack than works on all of the numbers
                // though I'm giving up some precision (anything beyond 9 digits is lost)
                //it keeps me from printing crap like 1.0000000000000000004 on the screen, though
                //and that's not real anyway because that's just the limitation 
                //of the hardware expressing decimal numbers as binary
            }
        }

        Matrix Result = new Matrix(result, A.getRows(), B.getColumns());

        return Result;
    }

    public Matrix mult(Matrix B) throws DimensionMismatchException {
        return mult(this, B);
    }

    public static Matrix zeroes(int rows, int columns) {
        double[][] zeroMatrix = new double[rows][columns];

        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < columns; ++j) {
                zeroMatrix[i][j] = 0;
            }
        }
        return new Matrix(zeroMatrix, rows, columns);
    }
}
