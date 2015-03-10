/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package matrix;

import java.util.logging.Level;
import java.util.logging.Logger;

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

        Matrix testMatrix4 = new Matrix(myArray, myArray.length, myArray[0].length);

        System.out.println("My testMatrix4 :");
        System.out.println(testMatrix4);

        double[][] mySecondArray = {
            {1, 2, 3},
            {4, 12, 6},
            {7, 8, 11}
        };

        Matrix testMatrix5 = new Matrix(mySecondArray, mySecondArray.length, mySecondArray[0].length);

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
            LUPDecomposition testDecompose = new LUPDecomposition(testSquare1);
        } catch (NonInvertibleMatrixException ex) {
            System.out.println("LUDecomposition of testDecompose object threw exception:\n" + ex.getMessage());
        }

        double[][] myThirdArray = {
            {7, 4, 2, 0},
            {6, 3, -1, 2},
            {4, 6, 2, 5},
            {8, 2, -7, 1}
        };

        SquareMatrix testSquare2 = new SquareMatrix(new Matrix(myThirdArray, 4, 4));

        try {
            LUPDecomposition testDecompose2 = new LUPDecomposition(testSquare2);

        } catch (NonInvertibleMatrixException ex) {
            System.out.println("LUDecomposition of testDecompose2 object threw exception:\n" + ex.getMessage());
        }

        System.out.println("Smoke test for testSquare2.determinant():\t" + testSquare2.determinant());
        System.out.println("Second smoke test for testSquare1.determinant():\t" + testSquare1.determinant());
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
        
        
    }

}
