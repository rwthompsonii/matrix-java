/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package matrix;

/**
 * class provides QR decomposition via householder reflections method
 * @author Ron
 */
public class QRDecomposition {
    public Matrix Q;
    public Matrix R;
    
    public QRDecomposition(Matrix A){
        Q = Matrix.zeroes(A.getRows(), A.getColumns());
        R = Matrix.zeroes(A.getRows(), A.getColumns());
    }
    
    public void decompose(Matrix A){
        Matrix q = new Matrix(Q);
    }
}
