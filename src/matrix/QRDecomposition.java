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

/**
 * class provides QR decomposition via householder reflections method
 *
 * @author Ron
 */
public class QRDecomposition implements MatrixConstants {

    public Matrix Q;
    public Matrix R;
    public int iterations;

    public QRDecomposition() {

    }

    public Matrix hessenberg(Matrix A) {

        int m = A.getRows();
        int n = A.getColumns();
        
        if (n <= 2) {
            //deals with pathological cases
            return A;
        }

        double[][] h = new double[m][n];
        
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                h[i][j] = A.getMatrix()[i][j];
            }
        }

        //first, we want to iterate over the columns of h from k = 0 to n-3
        //that's two columns from the end 
        for (int k = 0; k < n - 2; ++k) {
            System.out.println("iteration #" + (k+1));
            
            //move the kth column of h into v
            double[] v = new double[m - (k + 1)];
            //move the elements in the kth column
            //that I want zeroed into the v array
            for (int j = k + 1; j < m; ++j) {
                v[j - (k + 1)] = h[j][k];
            }
            //calculate the norm
            double alpha = -vector_norm(v);

            //hack for the sign function over the basis vector
            if (v[0] < 0) {
                alpha *= -1;
            }
            v[0] -= alpha;

            //scale the vector by its norm
            v = vector_scale(v, (1 / vector_norm(v)));
            
            
            Matrix V = new Matrix(v, ROW_VECTOR);
            Matrix Vt = new Matrix(v, COLUMN_VECTOR);

            //copy the submatrix into the temp matrix
            double[][] temp = new double[m - (k + 1)][n - (k + 1)];
            
            for (int i = k + 1; i < m; ++i) {
                for (int j = k + 1; j < n; ++j) {
                    temp[i - (k + 1)][j - (k + 1)] = h[i][j];
                }
            }

            Matrix B = new Matrix(temp);
            //temp = temp - 2 * Vt * (V * temp)
            try {
                B = B.subtract(Vt.scale(2).mult(V.mult(B)));
            } catch (DimensionMismatchException ex) {
                ex.printStackTrace(System.out);
                System.out.println("An unexpected exception occurred on first hessenberg line.  Bailing");
                System.exit(-1);
            }

            temp = B.getMatrix();

            //put it back in h, just reverse the assignment from before
            for (int i = k + 1; i < m; ++i) {
                for (int j = k + 1; j < n; ++j) {
                    h[i][j] = temp[i - (k + 1)][j - (k + 1)];
                }
            }

            //fixing roundoff errors from multiplication
            h[k + 1][k] = alpha;
            for (int i = k + 2; i < m; ++i) {
                h[i][k] = 0;
            }
            
            temp = new double[m][n - (k + 1)];

            //construct another temp matrix 
            for (int i = 0; i < m; ++i) {
                for (int j = k + 1; j < n; ++j) {
                    temp[i][j - (k + 1)] = h[i][j];
                }
            }

            B = new Matrix(temp);
            //Matrix C = minor(new Matrix(h), 0, k+1);
            
            
            
            try {
                B = B.subtract((B.mult(Vt)).mult(V).scale(2));
            } catch (DimensionMismatchException ex) {
                System.out.println("An unexpected exception occurred on second hessenberg lines.  Bailing");
                System.exit(-1);
            }
            
            temp = B.getMatrix();
            
            //put temp back into h and start the loop back over
            
            for (int i = 0; i < m; ++i) {
                for (int j = k + 1; j < n; ++j) {
                    h[i][j] = temp[i][j - (k + 1)] ;
                }
            }
            
        }

        return new Matrix(h);
    }

    public void decompose(Matrix A) {
        //initialize Q & R with sizes of A
        Q = Matrix.zeroes(A.getRows(), A.getColumns());
        R = Matrix.zeroes(A.getRows(), A.getColumns());

        Matrix[] qArray = new Matrix[A.getRows()];
        Matrix z = new Matrix(A);
        Matrix z1 = null;

        for (int k = 0; k < A.getRows() - 1 && k < A.getColumns(); ++k) {
            double[] e = new double[A.getRows()];
            double[] x = new double[A.getRows()];
            double a = 0.0;

            z1 = minor(z, k, k);
            z = z1;

            //move the kth column of z into x;
            for (int i = 0; i < z.getRows(); ++i) {
                x[i] = z.getMatrix()[i][k];
            }

            a = vector_norm(x);

            if (A.getMatrix()[k][k] > 0) {
                a *= -1;
            }

            for (int i = 0; i < A.getRows(); ++i) {
                if (k == i) {
                    e[i] = 1;
                } else {
                    e[i] = 0;
                }
            }

            //e = x + a * e
            for (int i = 0; i < A.getRows(); ++i) {
                e[i] = x[i] + a * e[i];
            }

            //e = e/norm(e)
            e = vector_scale(e, (1 / vector_norm(e)));

            //Q = I - 2*v*v^T where v=e
            Matrix Identity = Matrix.IdentityMatrix(A.getRows());
            Matrix Vector = new Matrix(vector_mult(e, e));
            try {
                qArray[k] = Identity.subtract(Vector.scale(2));
            } catch (DimensionMismatchException ex) {
                System.out.println("An unexpected exception occurred during I - 2*v*v^T, bailing.");
                System.exit(-1);
            }
            // z1 = q[k] * z
            try {
                z1 = qArray[k].mult(z);
            } catch (DimensionMismatchException ex) {
                System.out.println("An unexpected exception occurred during z1 = q[k]*z, bailing.");
                System.exit(-1);
            }
            z = z1;

        }

        Q = qArray[0];
        try {
            // R = Q * A
            R = Q.mult(A);
        } catch (DimensionMismatchException ex) {
            System.out.println("An unexpected exception occurred during R = Q * A, bailing.");
            System.exit(-1);
        }

        //Q^T = q[0]*q[1]*...*q[i] 
        for (int i = 1; i < A.getColumns() && i < A.getRows() - 1; ++i) {
            try {
                z1 = qArray[i].mult(Q);
            } catch (DimensionMismatchException ex) {
                System.out.println("An unexpected exception occurred during z1 = q[i]*Q, bailing.");
                System.exit(-1);
            }
            Q = z1;
        }

        try {
            //z = Q * A
            z = Q.mult(A);
        } catch (DimensionMismatchException ex) {
            System.out.println("An unexpected exception occurred during z = Q * A, bailing.");
            System.exit(-1);
        }
        //finval value of R is z
        R = z;
        //Q^T = final Q
        Q = Matrix.transpose(Q);

    }

    /**
     *
     * internal helper method used in decompose()
     *
     * @param A - the matrix to compute the sub-matrix of
     * @param row - the row that the sub-matrix will start on
     * @return - a Matrix populated with the identity matrix before A[row][row]
     * and data from the original matrix on/after
     */
    private Matrix minor(Matrix A, int row, int column) {
        double[][] result = new double[A.getRows()][A.getColumns()];

        for (int i = 0; i <= row; ++i) {
            result[i][i] = 1;
        }
        for (int i = row; i < A.getRows(); ++i) {
            for (int j = column; j < A.getColumns(); ++j) {
                result[i][j] = A.getMatrix()[i][j];
            }
        }

        return new Matrix(result);
    }

    /**
     *
     * internal helper method used in decompose()
     *
     * @param A - the matrix to compute the sub-matrix of
     * @param row - the row that the sub-matrix will start on
     * @return - a Matrix populated with the identity matrix before A[row][row]
     * and data from the original matrix on/after
     */
    private double vector_norm(double[] vector) {
        double sum = 0;
        int n = vector.length;
        for (int i = 0; i < n; ++i) {
            sum += vector[i] * vector[i];
        }

        return Math.sqrt(sum);
    }

    /**
     *
     * internal helper method used in decompose()
     *
     * @param vector - the vector that will be scaled - note it is modified
     * @param scale - the scaling factor to be applied where v = v * s
     * @return - the original vector that has been scaled
     */
    private double[] vector_scale(double[] vector, double scale) {
        for (int i = 0; i < vector.length; ++i) {
            vector[i] *= scale;
        }
        return vector;
    }

    private double[][] vector_mult(double[] a, double[] b) {
        if (a.length != b.length) {
            throw new IllegalArgumentException("vector multiplication can only be applied on vectors of equal length");
        }

        double[][] result = new double[a.length][a.length];

        for (int i = 0; i < a.length; ++i) {
            for (int j = 0; j < a.length; ++j) {
                result[i][j] = a[i] * b[j];
            }
        }

        return result;
    }
}
