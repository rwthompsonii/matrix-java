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

    //note that equals does NOT compare with equality operator, it uses EPSILON for "close enough"
    public boolean equals(Matrix other) {
        if (other.getClass() != getClass() || other.getColumns() != getColumns() || other.getRows() != getRows()) {
            return false;
        } else {
            double[][] otherMatrix = other.getMatrix();
            for (int i = 0; i < rows; ++i) {
                for (int j = 0; j < columns; ++j) {
                    if (Math.abs(otherMatrix[i][j] - matrix[i][j]) < EPSILON) {
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
            }
        }

        Matrix Result = new Matrix(result, A.getRows(), B.getColumns());

        return Result;
    }
    
    public static Matrix transpose(Matrix A) {
        double[][] result = new double[A.getColumns()][A.getRows()];
        
        double[][] a = A.getMatrix();
        
        for(int i = 0; i < A.getRows(); ++i){
            for(int j = 0; j < A.getColumns(); ++j){
                result[i][j] = a[j][i];
            }
        }
  
        Matrix Result = new Matrix(result, A.getColumns(), A.getRows());
        
        return Result;
    }
    
    public Matrix transpose() {
        return transpose(this);
    }
    
    public static Matrix scale(Matrix A, double v) {
    
        double[][] a = A.getMatrix().clone();
        
        for(int i = 0; i < A.getRows(); ++i) {
            for(int j = 0; j < A.getColumns(); ++j) {
                a[i][j] *= v;
            }
        }
        
        Matrix Result = new Matrix(a, A.getRows(), A.getColumns());
    
        return Result;
    }
    
    public Matrix scale(double v) {
        return scale(this, v);
    }
    
    public static Matrix add(Matrix A, Matrix B) throws DimensionMismatchException {
        if(A.getRows() != B.getRows() || A.getColumns() != B.getColumns()) {
            throw new DimensionMismatchException("This operation can only be performed on matrices of equal dimensions.");
        }
        
        double[][] result = A.getMatrix().clone();
        
        for(int i = 0; i < A.getRows(); ++i) {
            for(int j = 0; j < A.getColumns(); ++j) {
                result[i][j] = A.getMatrix()[i][j] + B.getMatrix()[i][j];
            }
        }
        
        Matrix Result = new Matrix(result, A.getRows(), A.getColumns());
        
        return Result;
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
