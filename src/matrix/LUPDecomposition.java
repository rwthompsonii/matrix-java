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
public class LUPDecomposition {

    public SquareMatrix U;
    public SquareMatrix L;
    public SquareMatrix P;
    public int permutations;

    public LUPDecomposition() {
    }

    /**
     *
     * @param A - the input Matrix
     * @return returns number of permutations, useful for calculating
     * determinants. ignore it if you don't need it
     * @throws matrix.NonInvertibleMatrixException
     */
    public int decompose(SquareMatrix A) throws NonInvertibleMatrixException {
        U = new SquareMatrix(Matrix.zeroes(A.getRows(), A.getColumns()));
        L = new SquareMatrix(Matrix.IdentityMatrix(A.getRows()));
        SquareMatrix aPrime = null;
        pivot(A);
        try {
            aPrime = new SquareMatrix(P.mult(A));
        } catch (DimensionMismatchException ex) {
            System.out.println("Exception: " + ex.getMessage());
            System.exit(-1);//this really shouldn't happen so i'm going to bail here
        }

        for (int i = 0; i < A.getRows(); ++i) {
            for (int j = 0; j < A.getColumns(); ++j) {
                double s = 0.0;
                if (j <= i) {
                    for (int k = 0; k < j; ++k) {
                        s += L.getMatrix()[j][k] * U.getMatrix()[k][i];
                    }
                    U.getMatrix()[j][i] = aPrime.getMatrix()[j][i] - s;
                    if (i == j) {
                        if (Math.abs(U.getMatrix()[j][i]) <= 0.000001) {
                            L = null;
                            U = null;
                            P = null;//cleaning these up so they don't get used
                            throw new NonInvertibleMatrixException("Matrix cannot be decomposed, it is singular.");
                        }
                    }
                }

                s = 0.0;
                if (j > i) {
                    for (int k = 0; k < i; ++k) {
                        s += L.getMatrix()[j][k] * U.getMatrix()[k][i];
                    }
                    L.getMatrix()[j][i] = (aPrime.getMatrix()[j][i] - s) / U.getMatrix()[i][i];
                }
            }
        }
        /*
         System.out.println("\nA:\n" + A + "\naPrime:\n" + aPrime + "\nL:\n" + L + "\nU:\n" + U + "\nP:\n" + P);
         try {
         System.out.println("Result of P.mult(A):\n" + P.mult(A));
         System.out.println("Result of L.mult(U):\n" + L.mult(U));
         } catch (DimensionMismatchException ex) {
         ex.printStackTrace();
         System.exit(-1);//this really shouldn't happen so i'm going to bail here
         }*/

        return permutations;
    }
    
    /**
     *
     * @param a - the matrix to be decomposed, note that if it is not square an error will occur, calling code must verify this.  
     * if 'a' is not invertible an exception will be thrown
     * @param l - output l matrix, lower unit triangular
     * @param u - output u matrix, upper triangular
     * @param p - output permutation matrix
     * @return - number of row permutations in p matrix from identity matrix, useful for calculating determinants
     * @throws NonInvertibleMatrixException
     */
    public static int decompose (double[][] a, double[][] l, double[][] u, double[][] p) throws NonInvertibleMatrixException {
        LUPDecomposition LUP = new LUPDecomposition();
        SquareMatrix A = new SquareMatrix(a);
        
        int result =  LUP.decompose(A);
        
        l = LUP.L.getMatrix();
        u = LUP.U.getMatrix();
        p = LUP.P.getMatrix();
        
        return result;
    }

    //helper function for decompose

    private void pivot(SquareMatrix A) {
        permutations = 0;
        P = new SquareMatrix(Matrix.IdentityMatrix(A.getRows()));

        for (int i = 0; i < A.getRows(); ++i) {
            int maxJ = i;

            for (int j = i; j < A.getRows(); ++j) {
                if (Math.abs(A.getMatrix()[j][i]) > Math.abs(A.getMatrix()[maxJ][i])) {
                    maxJ = j;
                }
            }

            if (maxJ != i) {
                double temp;
                for (int k = 0; k < A.getRows(); ++k) {
                    temp = P.getMatrix()[i][k];
                    P.getMatrix()[i][k] = P.getMatrix()[maxJ][k];
                    P.getMatrix()[maxJ][k] = temp;
                }
                permutations++;
            }
        }
    }
}
