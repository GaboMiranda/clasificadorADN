/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.mycompany.biogabo;

import java.io.IOException;
import java.util.ArrayList;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;

/**
 *
 * @author morelos
 */
public class ClasificadorBayes {
    FastaIO io =new FastaIO();
    
    static String []combinaciones=new String[125]; 
    static private String tiplete;
    static ArrayList<Integer> valores = new ArrayList();
    public void crearTripletes() throws IOException, CompoundNotFoundException
    {
        String bases="ATGCN";
        int cont=0;
        for (int i = 0; i < bases.length(); i++) 
        {
            for (int j = 0; j < bases.length(); j++) 
            {
                for (int k = 0; k < bases.length(); k++) 
                {
                    combinaciones[cont]=bases.charAt(i)+""+bases.charAt(j)+""+bases.charAt(k);
                    cont++;
                }                
            }
            
        }
    }
    public int compararTrip(String seq) throws IOException, CompoundNotFoundException
    {
        this.crearTripletes();
        
        seq=seq.substring(0, 15);
        for (int i = 0; i < seq.length(); i++) 
        {
            String corte=seq.substring(i,i+3);
            for (int j = 0; j < combinaciones.length; j++) 
            {
                if(corte.compareTo(combinaciones[j])==0)
                {
                    valores.add(j);
//                    System.out.print(" "+j);
                }
            }
            i+=2;
        }
        return this.evaluarTrip();
    }
    
    public int evaluarTrip() throws IOException{
        ArrayList<Double> probS = new ArrayList();
        ArrayList<Double> resultados = new ArrayList();
        double prob = 1.0;
        double probG = 1.0;
        int mGrupo=0;
        FastaIO io = new FastaIO();
        int ngrupos=0;
        double dngrupos=io.nGrupos(0);
        ngrupos=(int)dngrupos;
        for (int i = 1; i <= ngrupos; i++) 
        {
            for (int j = 0; j < valores.size(); j++) 
            {
//                System.out.println();
                
                prob*=io.ReadProbAtrb(i, valores.get(j),j+1);
//                System.out.println(prob);
//                System.out.println(prob+" "+io.ReadProbAtrb(i, valores.get(j)));
            }
            probG=io.nGrupos(i);
//            System.out.println(prob+" "+probG);
            resultados.add(prob*probG);
            prob=1;
        }
        mGrupo=this.mayorProb(resultados);
        System.out.println("La secuencia pertenece al Grupo: "+mGrupo+" "+resultados.size());
        valores.clear();
        return mGrupo;
    }
    public int mayorProb(ArrayList<Double> resultados)
    {
        double iNumeroMayor=resultados.get(0);
        int iPosicion = 0;
        for (int i = 0; i <resultados.size(); i++) 
        {
//            System.out.println(resultados.get(i));
            if (resultados.get(i)>iNumeroMayor)
            {
                iNumeroMayor = resultados.get(i);
                iPosicion = i;
            } 
        }
        return iPosicion+1;
    }
    
}
