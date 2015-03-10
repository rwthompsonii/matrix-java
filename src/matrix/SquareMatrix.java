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
public class SquareMatrix extends Matrix {

    private static LUPDecomposition LUP;

    public SquareMatrix() {
        super();
        assert (rows == columns);//this one is only necessary if 
        //i change the default ctor on super and forget to change this
        //, so it's more of a reminder
    }

    public SquareMatrix(int[][] matrix, int size) {
        super(matrix, size, size);
        assert (matrix[0].length == matrix.length);
        assert (matrix.length == size);
    }

    public SquareMatrix(float[][] matrix, int size) {
        super(matrix, size, size);
        assert (matrix[0].length == matrix.length);
        assert (matrix.length == size);
    }

    public SquareMatrix(double[][] matrix, int size) {
        super(matrix, size, size);
        assert (matrix[0].length == matrix.length);
        assert (matrix.length == size);
    }

    public SquareMatrix(SquareMatrix square) {
        super(square.getMatrix(), square.getRows(), square.getColumns());
        //no check since it's copying an existing square
    }

    public SquareMatrix(Matrix maybe) {
        this(maybe.getMatrix(), maybe.getRows());
    }

    /**
     *
     * @param A a square matrix to calculate the determinant of, which is
     * returned
     * @return
     */
    public static double determinant(SquareMatrix A) {
        try {
            LUP = new LUPDecomposition(A);
        } catch (NonInvertibleMatrixException ex) {
            LUP = null;
            return 0;
        }
        double det = 1;
        double[][] u = LUP.U.getMatrix();
        /*
         System.out.println("\nA:\n" + A + "\nL:\n" + L + "\nU:\n" + U + "\nP:\n" + P + "\nperms:\n" + perms);
         try {
         System.out.println("Result of P.mult(A):\n" + P.mult(A));
         System.out.println("Result of L.mult(U):\n" + L.mult(U));
         } catch (DimensionMismatchException ex) {
         ex.printStackTrace();
         exit(-1);//this really shouldn't happen so i'm going to bail here
         }*/

        //still need to account for P matrix.  --done
        //L matrix is not used because it is a unit lower triangular matrix, 
        //and its main diagonal entries are all 1
        for (int i = 0; i < LUP.U.getRows(); ++i) {
            det *= u[i][i];
        }

        if (LUP.permutations % 2 != 0) {
            det *= -1;
        }

        return det;
    }

    /**
     *
     * @return Determinant of matrix in class
     */
    public double determinant() {
        return determinant(this);
    }

    public static SquareMatrix invert(SquareMatrix A) throws NonInvertibleMatrixException {
        double detA = determinant(A);
        if (0 == detA) {
            throw new NonInvertibleMatrixException("Cannot invert matrix, determinant is 0.");
        }
        int n = A.getRows();
        double[][] inverse = new double[n][n];
        
        for(int i = 0; i < n; ++i) {
            double [] columnSolved = solve(i);
            for (int j = 0; j < n; ++j) {
                inverse[j][i] = columnSolved[j];
            }
        }

        SquareMatrix Inverse = new SquareMatrix(inverse, n);

        return Inverse;
    }

    //internal helper method used in invert()
    private static double[] solve(int column) {
        int n = LUP.L.getRows();

        double[] result = new double[n];

        for (int i = 0; i < n; ++i) {
            result[i] = LUP.P.getMatrix()[i][column];
        }
        
        for (int i = 1; i < n; ++i) {
            double s = result[i];
            for(int j = 0; j < i; ++j) {
                s -= result[j] * LUP.L.getMatrix()[i][j];
            }
            result[i] = s;
        }
        
        result[n-1] /= LUP.U.getMatrix()[n-1][n-1];
        
        for (int i = n-2; i >= 0; --i) {
            double s = result[i];
            for (int j = i +1; j < n; ++j) {
                s -= result[j] * LUP.U.getMatrix()[i][j];
            }
            result[i] = s/LUP.U.getMatrix()[i][i];
        }

        return result;
    }
}
