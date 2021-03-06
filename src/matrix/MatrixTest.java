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

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.util.Precision;

/**
 *
 * @author ron
 */
public class MatrixTest {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        Matrix testMatrix = new Matrix();
        System.out.println("My testMatrix (should be identity 2x2):");
        System.out.println(testMatrix);

        Matrix testMatrix2 = Matrix.IdentityMatrix(3);
        System.out.println("My testMatrix2 (should be identity 3x3) : ");
        System.out.println(testMatrix2);

        Matrix testMatrix3 = Matrix.IdentityMatrix(3);
        System.out.println("My testMatrix3 (should be identity 3x3) : ");
        System.out.println(testMatrix3);

        int[][] myArray = {
            {1, 0, 0},
            {0, 1, 0},
            {0, 0, 2}
        };

        Matrix testMatrix4 = new Matrix(myArray);

        System.out.println("My testMatrix4 :");
        System.out.println(testMatrix4);

        double[][] mySecondArray = {
            {1, 2, 3},
            {4, 12, 6},
            {7, 8, 11}
        };

        Matrix testMatrix5 = new Matrix(mySecondArray);

        System.out.println("My testMatrix5 :");
        System.out.println(testMatrix5);

        System.out.println("Result of testMatrix2.equals(testMatrix3): " + testMatrix2.equals(testMatrix3));
        System.out.println("Result of testMatrix.equals(testMatrix3): " + testMatrix.equals(testMatrix3));
        System.out.println("Result of testMatrix3.equals(testMatrix4): " + testMatrix3.equals(testMatrix4));
        try {
            System.out.println("\nResult of Matrix.mult(testMatrix4, testMatrix3):\n" + Matrix.mult(testMatrix4, testMatrix3));
            System.out.println("\nResult of testMatrix3.mult(testMatrix4):\n" + testMatrix3.mult(testMatrix4));
            System.out.println("\nResult of testMatrix4.mult(testMatrix5):\n" + testMatrix4.mult(testMatrix5));
            System.out.println("\nResult of testMatrix5.mult(testMatrix4):\n" + testMatrix5.mult(testMatrix4));
        } catch (DimensionMismatchException ex) {
            System.out.println(ex.getMessage());
        }

        SquareMatrix testSquare1 = new SquareMatrix(testMatrix5);
        try {
            LUPDecomposition testDecompose = new LUPDecomposition();
            testDecompose.decompose(testSquare1);
        } catch (NonInvertibleMatrixException ex) {
            System.out.println("LUDecomposition of testDecompose object threw exception:\n" + ex.getMessage());
        }

        double[][] myThirdArray = {
            {7, 4, 2, 0},
            {6, 3, -1, 2},
            {4, 6, 2, 5},
            {8, 2, -7, 1}
        };

        SquareMatrix testSquare2 = new SquareMatrix(new Matrix(myThirdArray));

        try {
            LUPDecomposition testDecompose2 = new LUPDecomposition();
            testDecompose2.decompose(testSquare2);
        } catch (NonInvertibleMatrixException ex) {
            System.out.println("LUDecomposition of testDecompose2 object threw exception:\n" + ex.getMessage());
        }

        System.out.println("Smoke test for testSquare2.determinant():\t" + testSquare2.determinant());
        System.out.println("Second smoke test for testSquare1.determinant():\t" + testSquare1.determinant() + "\n");
        try {
            SquareMatrix InvertedTestSquare2 = SquareMatrix.invert(testSquare2);
            System.out.println("InvertedTestSquare2:\n" + InvertedTestSquare2);
            System.out.println("Result of testSquare2.mult(InvertedTestSquare2): (should be identity matrix)\n" + testSquare2.mult(InvertedTestSquare2));
            System.out.println("Result of InvertedTestSquare2.mult(TestSquare2): (should be identity matrix)\n" + InvertedTestSquare2.mult(testSquare2));
        } catch (NonInvertibleMatrixException ex) {
            System.out.println("Inverting matrix testSquare2 threw exception:\n" + ex.getMessage());
        } catch (DimensionMismatchException ex) {
            System.out.println("Multiplying Matrix failed. (shouldn't have happened, exiting");
            System.exit(-1);
        }

        SquareMatrix testSquare3 = new SquareMatrix(Matrix.transpose(testSquare2));

        System.out.println("testSquare2:\n" + testSquare2 + "\ntestSquare3: (transpose of testSquare2)\n" + testSquare3);
        try {
            System.out.println("Result of testSquare2.subtract(testSquare3)\n" + testSquare2.subtract(testSquare3));
        } catch (DimensionMismatchException ex) {
            System.out.println("Subtracting Matrix failed. (shouldn't have happened, exiting");
            System.exit(-1);
        }

        /*float[][] fail = {
         {1, 2, 3},
         {4, 5, 7},
         {3, 6, 8},
         {4, 5, 9}
         };*/
        //Matrix goodMatrix = new Matrix(fail);
        //SquareMatrix badMatrix = new SquareMatrix(goodMatrix);
        //System.out.print(badMatrix);
        Matrix random = Matrix.random(3, 5, 100, 1);

        System.out.println("random matrix: (should be 3x5)\n" + random);

        int BMFsize = 500;
        System.out.println("smoke checking my mult method with " + BMFsize + "x" + BMFsize + " matrices.\n");

        Matrix big1 = Matrix.random(BMFsize, BMFsize, 2, -2);
        Matrix big2 = Matrix.random(BMFsize, BMFsize, 2, -2);
        SquareMatrix big3 = null;

        try {
            big3 = new SquareMatrix(big1.mult(big2));
        } catch (DimensionMismatchException ex) {
            System.out.println("this exception shouldn't have happened, bailing");
            System.exit(-1);
        }

        System.out.println("smoke checking determinant method with " + BMFsize + "x" + BMFsize + " matrix.\n");
        System.out.println("Result of big3.determinant():\t" + big3.determinant() + "\n");

        System.out.println("Performing basic check of QR Decomp with known 4x4 matrix.\n");

        double[][] qrTest = {
            {4, 8, 9, 1},
            {2, -3, 1, -1},
            {8, 15, 3, 1},
            {7, 1, 2, 3}};

        Matrix QRTest = new Matrix(qrTest);

        QRDecomposition qr = new QRDecomposition();

        qr.decompose(QRTest);

        System.out.println("\nQ:\n" + qr.Q + "\nR:\n" + qr.R + "\nQRtest: (original matrix)\n" + QRTest);

        try {
            System.out.println("Result of Q.mult(R): (should be the original matrix)\n" + qr.Q.mult(qr.R));
        } catch (DimensionMismatchException ex) {
            System.out.println("this exception shouldn't have happened, bailing");
            System.exit(-1);
        }

        try {
            System.out.println("Result of Q.mult(Q.transpose()): (should be the identity matrix)\n" + qr.Q.mult(qr.Q.transpose()));
        } catch (DimensionMismatchException ex) {
            System.out.println("this exception shouldn't have happened, bailing");
            System.exit(-1);
        }

        /*System.out.println("Beginning initial testing of eigenvalue function on QRTest:\n" + QRTest);

         Complex[] eigen = SquareMatrix.eigenvalues(new SquareMatrix(QRTest));

         for (Complex e : eigen) {
         System.out.println("eigenvalue:" + e);
         }*/
        int SMFSize = 10;

        //SquareMatrix big4 = new SquareMatrix(Matrix.random(SMFSize, SMFSize, -5, 5));
        //this matrix is causing a bug, using it to troubleshooot
        double[][] big4Matrix = {
            {-2.476, -2.814, 4.29, -3.649},
            {2.839, -2.859, 1.623, -2.926},
            {-0.392, -3.206, -0.401, -2.174},
            {2.241, -4.435, -3.963, 4.102}};
        SquareMatrix big4 = new SquareMatrix(big4Matrix);

        if (SMFSize < 6) {
            System.out.println("Smoke test of eigenvalue function with " + SMFSize + "x" + SMFSize + " matrix:\n" + big4);
        }//only print it if it will fit on the screen
        Complex eigen2[] = big4.eigenvalues();

        //System.out.println(big4.toCopyableString());//used for debugging
        int eig = 1;
        for (Complex e : eigen2) {
            System.out.println("eigenvalue #" + eig++ + ":\t" + Precision.round(e.getReal(), 3) + " + " + Precision.round(e.getImaginary(), 3) + "i");
        }

        double[] vector = new double[5];
        double[] vector2 = new double[5];

        Matrix v1 = new Matrix(vector, MatrixConstants.ROW_VECTOR);
        Matrix v2 = new Matrix(vector2, MatrixConstants.COLUMN_VECTOR);

        System.out.println("\nv1: (should be 1x5)\n" + v1 + "\nv2: (should be 5x1)\n" + v2);

        double[][] qrTest2 = {
            {4, 2, 2, 1},
            {2, -3, 1, 1},
            {2, 1, 3, 1},
            {1, 1, 1, 2}};

        Matrix QRTest2 = new Matrix(qrTest2);

        System.out.println("Testing hessenberg() with QRTest2:\n" + QRTest2);

        System.out.println("Hessenberg output for QRTest2:\n" + qr.hessenberg(QRTest2));

        int big6size = 5;
        SquareMatrix big6 = new SquareMatrix(Matrix.random(big6size, big6size, -5, 5));

        System.out.println("Testing hessenberg() with big6:\n" + big6);

        System.out.println("Hessenberg output with big6:\n" + qr.hessenberg(big6));

        //System.out.println(big6.toCopyableString());
    }

}
