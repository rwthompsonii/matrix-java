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
class DimensionMismatchException extends Exception {

    public DimensionMismatchException() {
        super();
    }
    public DimensionMismatchException(String message){
        super(message);
    }
    
}
