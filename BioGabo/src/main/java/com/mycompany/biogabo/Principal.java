/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.mycompany.biogabo;
import java.io.IOException;
import java.util.Scanner;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
/**
 *
 * @author morelos
 */
public class Principal 
{
    public  static void main(String[] args) throws IOException, CompoundNotFoundException 
    {
        interfaz in = new interfaz();
        in.setVisible(true);
        Scanner teclado= new Scanner(System.in);
        int opc;
        do{
            System.out.println("\nElegir una opción:");
            System.out.println("1.-Entrenar algoritmo");
            System.out.println("2.-Clasificar secuencia");
            System.out.println("3.-Salir");
            System.out.print("R = ");
            opc = teclado.nextInt();
            switch(opc)
            {
                case 1:
                    
                    break;
                case 2:
                    
                    Scanner te= new Scanner(System.in);
                    System.out.println("Ingresa la secuencia a agrupar:");
                    
                    
                    break;
                case 3:
                    break;
                default: 
                    System.out.println("opción no valida");
                    break;
            }
        }while(opc!=3);
    }
}
