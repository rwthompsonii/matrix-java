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
