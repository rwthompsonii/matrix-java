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

import java.math.BigDecimal;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.math3.complex.Complex;

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

    public Complex[] eigenvalues() {
        return eigenvalues(this);
    }

    public static Complex[] eigenvalues(SquareMatrix A) {
        Complex[] e = new Complex[A.getRows()];

        QRDecomposition qr = new QRDecomposition();

        qr.iterations = 0;
        int total_iter = 0;
        int num_eigen_found = 0;

        //in general, QR decomposition will converge faster from an upper
        //Hessenberg matrix.  so, first things first, we bring QRIterator to that form
        SquareMatrix QRIterator = new SquareMatrix(qr.hessenberg(A));
        //SquareMatrix QRIterator = new SquareMatrix(A);

        int max = MAX_ITERATIONS;
        double lastElement;
        SquareMatrix ScaledIdentity;
        do {

            System.out.println("Pre-decompose: QRIterator (Iteration#" + (qr.iterations + 1) + "):\n" + QRIterator);
            if (QRIterator.getRows() == 1) {
                //very last 1x1 element in matrix
                e[num_eigen_found++] = new Complex(
                        QRIterator.getMatrix()[0][0]
                );
                break;
            } else {

                lastElement = QRIterator.getMatrix()[QRIterator.getRows() - 1][QRIterator.getColumns() - 1];
                ScaledIdentity = new SquareMatrix(Matrix.IdentityMatrix(QRIterator.getRows()).scale(lastElement));
                try {
                    QRIterator = new SquareMatrix(QRIterator.subtract(ScaledIdentity));
                } catch (DimensionMismatchException ex) {
                    System.out.println("Unexpected execption during QRIterator -= I*alpha, bailing.");
                    System.exit(-1);

                }
                qr.decompose(QRIterator);
            }
            try {
                QRIterator = new SquareMatrix(qr.R.mult(qr.Q)/*.add(ScaledIdentity)*/);

            } catch (DimensionMismatchException ex) {
                System.out.println("An unexpected exception occurred during QRIterator = R*Q, bailing.");
                System.exit(-1);
            }
            qr.iterations++;

            //testing indicates that MAX_ITERATIONS iterations should be more than sufficient to converge, if its going to at all
            if (qr.iterations == max || Math.abs(QRIterator.getMatrix()[QRIterator.getRows() - 1][QRIterator.getColumns() - 2]) < CONVERGENCE_CHECK) {
                System.out.println("QRIterator (at max iteration or converged) (Iteration#" + qr.iterations + "):\n" + QRIterator + "\nlastElement value:\t" + lastElement);
                if (Math.abs(QRIterator.getMatrix()[QRIterator.getRows() - 1][QRIterator.getColumns() - 2]) < CONVERGENCE_CHECK) {
                    //then the value at M[n][n] is an eigenvalue and it is real
                    e[num_eigen_found++] = new Complex(
                            QRIterator.getMatrix()[QRIterator.getRows() - 1][QRIterator.getColumns() - 1]
                    );

                    //System.out.println("e[" + (num_eigen_found - 1) + "]:\t" + e[num_eigen_found - 1] + "\nQRIterator before deflation:\n" + QRIterator);
                    double[][] deflatedMatrix = deflate(QRIterator.getMatrix(), 1);
                    QRIterator = new SquareMatrix(deflatedMatrix);

                    total_iter += qr.iterations;
                    qr.iterations = 0;  //reset the iterations counter to find the next eigenvalue

                    //System.out.println("\nQRIterator after deflation:\n" + QRIterator);
                    /*if (2 <= QRIterator.getRows()) {
                     total_iter += qr.iterations;
                     qr.iterations = 0;  //reset the iterations counter to find the next eigenvalue
                     } else {
                     if (QRIterator.getRows() == 2) {
                     //i'm on the last 2x2 set and it's converged, so i'm going to pull the eigenvalues off the diagonal
                     e[num_eigen_found++] = new Complex(
                     QRIterator.getMatrix()[QRIterator.getRows() - 1][QRIterator.getColumns() - 1]
                     );
                     e[num_eigen_found++] = new Complex(
                     QRIterator.getMatrix()[QRIterator.getRows() - 2][QRIterator.getColumns() - 2]
                     );

                     } else if (QRIterator.getRows() == 1) {
                     //i'm on the very last 1x1 submatrix and it contains the last eigenvalue
                     e[num_eigen_found++] = new Complex(
                     QRIterator.getMatrix()[QRIterator.getRows() - 1][QRIterator.getColumns() - 1]
                     );
                     }

                     break;
                     }*/
                } else {
                    //this is a 2x2 matrix with either real or complex roots.  need to find them.
                    //characteristic equation of 2x2 array => E^2 - (w + z)E + (wz - xy) = 0 where E = eigenvalue (possibly pair, possibly singular, possibly real, possibly complex)
                    // and the matrix {{w, x}, {y, z}} is the input array, the task is to calculate the root(s) of that equation
                    //that is a quadratic equation => (root = (-b +- sqrt(b^2  - 4ac))/2a)
                    //determinant b^2 - 4ac will determine behavior of roots => positive means 2 real roots, 0 means 1 repeated real root, negative means conjugate pair of imaginary roots

                    //first, get the wxyz from the (possibly bigger) matrix
                    int n = QRIterator.getRows();
                    double w = QRIterator.getMatrix()[n - 2][n - 2];
                    double x = QRIterator.getMatrix()[n - 2][n - 1];
                    double y = QRIterator.getMatrix()[n - 1][n - 2];
                    double z = QRIterator.getMatrix()[n - 1][n - 1];

                    //a not used since it's = 1
                    double b = -(w + z);
                    double c = (w * z - x * y);

                    //calculate determinant of quadratic equation
                    double determ = b * b - 4 * c;

                    if (determ >= 0) {
                        //one or two real roots 
                        double sqrt_determ_real = Math.sqrt(determ);
                        e[num_eigen_found++] = new Complex((-b + sqrt_determ_real) / 2.0);
                        e[num_eigen_found++] = new Complex((-b - sqrt_determ_real) / 2.0);
                        //in the zero determinant case that's simply going to add the same eigenvalue to the list twice.  I'm ok with that for now.
                    } else if (determ < 0) {
                        //conjugate pair of complex roots
                        double sqrt_determ_imag = Math.sqrt(-determ);
                        e[num_eigen_found++] = new Complex(-b / 2.0, sqrt_determ_imag / 2.0);
                        e[num_eigen_found++] = new Complex(-b / 2.0, -sqrt_determ_imag / 2.0);
                    }

                    if (QRIterator.getRows() > 2) {
                        total_iter += qr.iterations;
                        qr.iterations = 0;  //reset the iterations counter to find the next eigenvalue
                        double[][] deflatedMatrix = deflate(QRIterator.getMatrix(), 2);
                        QRIterator = new SquareMatrix(deflatedMatrix);
                    }
                }
            }
            //QRIterator = new SquareMatrix(qr.hessenberg(QRIterator));

        } while (qr.iterations < max);

        //used for debugging here
        /*System.out.println("Finished iterating.  Iterations:\t" + qr.iterations
         + "\nFinal value of qr.Q:\n" + qr.Q + "\nFinal value of qr.R:\n" + qr.R
         + "\nFinal value of QRIterator:\n" + QRIterator
         + "\nOriginal SquareMatrix A:\n" + A);
         */
        return e;
    } //internal helper method called from eigenvalues function 

    private static boolean hasConverged(SquareMatrix A) {
        double[][] a = A.getMatrix().clone();//do not want to modify it here

        for (int i = 0; i < a.length; ++i) {
            for (int j = i - 1; j >= 0; --j) {
                if (MatrixConstants.CONVERGENCE_CHECK < Math.abs(a[i][j])) {
                    return false;
                }
            }
        }
        return true;
    }

    private static double[][] deflate(double[][] matrix, int numToDelete) {
        double[][] deflated = new double[matrix.length - numToDelete][matrix[0].length - numToDelete];

        for (int i = 0; i < deflated.length; ++i) {
            for (int j = 0; j < deflated[0].length; ++j) {
                deflated[i][j] = matrix[i][j];
            }
        }
        return deflated;
    }

}
