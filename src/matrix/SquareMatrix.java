/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package matrix;

import java.math.BigDecimal;

/**
 *
 * @author ron
 */
public class SquareMatrix extends Matrix {

    private static LUPDecomposition LUP;

    //I wanted to use asserts instead of exceptions but I found out through testing
    //that assert(expression) is not as meaninful to the java compiler/runtime as 
    //it is in c/c++ => i'm using an unchecked exception, because I *really*
    //don't want client code to continue execution after that exception is thrown

    public SquareMatrix() {
        super();
        if (rows != columns) {
            throw new IllegalArgumentException("Input matrix must have equal rows and columns.");
        }//this one is only necessary if 
        //i change the default ctor on super and forget to change this here
        //, so it's more of a reminder
    }

    public SquareMatrix(int[][] matrix) {
        super(matrix);
        if (rows != columns) {
            throw new IllegalArgumentException("Input matrix must have equal rows and columns.");
        }
    }

    public SquareMatrix(float[][] matrix) {
        super(matrix);
        if (rows != columns) {
            throw new IllegalArgumentException("Input matrix must have equal rows and columns.");
        }
    }

    public SquareMatrix(double[][] matrix) {
        super(matrix);
        if (rows != columns) {
            throw new IllegalArgumentException("Input matrix must have equal rows and columns.");
        }
    }

    public SquareMatrix(SquareMatrix square) {
        super(square);
        //no check since it's copying an existing square
    }

    public SquareMatrix(Matrix maybe) {
        this(maybe.getMatrix());
        //no check since it's going to be checked in the other ctor
    }

    /**
     *
     * @param A a square matrix to calculate the determinant of, which is
     * returned
     * @return
     */
    public static BigDecimal determinant(SquareMatrix A) {
        try {
            LUP = new LUPDecomposition();
            LUP.decompose(A);
        } catch (NonInvertibleMatrixException ex) {
            LUP = null;
            return new BigDecimal(0.0);
        }
        BigDecimal det = new BigDecimal(1);
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
            det = det.multiply(new BigDecimal(u[i][i]));
        }

        if (LUP.permutations % 2 != 0) {
            det = det.multiply(new BigDecimal(-1));
        }

        return det;
    }

    /**
     *
     * @return Determinant of matrix in class
     */
    public BigDecimal determinant() {
        return determinant(this);
    }

    public static SquareMatrix invert(SquareMatrix A) throws NonInvertibleMatrixException {
        BigDecimal detA = determinant(A);
        if (detA == new BigDecimal(0)) {
            throw new NonInvertibleMatrixException("Cannot invert matrix, determinant is 0.");
        }
        int n = A.getRows();
        double[][] inverse = new double[n][n];

        for (int i = 0; i < n; ++i) {
            double[] columnSolved = solve(i);
            for (int j = 0; j < n; ++j) {
                inverse[j][i] = columnSolved[j];
            }
        }

        return new SquareMatrix(inverse);
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
            for (int j = 0; j < i; ++j) {
                s -= result[j] * LUP.L.getMatrix()[i][j];
            }
            result[i] = s;
        }

        result[n - 1] /= LUP.U.getMatrix()[n - 1][n - 1];

        for (int i = n - 2; i >= 0; --i) {
            double s = result[i];
            for (int j = i + 1; j < n; ++j) {
                s -= result[j] * LUP.U.getMatrix()[i][j];
            }
            result[i] = s / LUP.U.getMatrix()[i][i];
        }

        return result;
    }
}
