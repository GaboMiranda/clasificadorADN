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
    }
}
