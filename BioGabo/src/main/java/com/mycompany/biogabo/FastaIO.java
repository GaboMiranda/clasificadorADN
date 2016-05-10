/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.mycompany.biogabo;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedHashMap;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.io.FastaReaderHelper;
/**
 *
 * @author morelos
 */
public class FastaIO {
    private final String url="/var/www/html/Gabo/";
    private LinkedHashMap<String, DNASequence> a;

    public LinkedHashMap<String, DNASequence> readFasta(String ruta,int tipo) throws IOException, CompoundNotFoundException {                 
        
    try {
        //Try with the FastaReaderHelper
        if(tipo==1)
        {
            File file = new File(url+ruta+".fasta");
            a = FastaReaderHelper.readFastaDNASequence(file);
        }
        else
        {
            File file = new File(ruta+".fasta");
            a = FastaReaderHelper.readFastaDNASequence(file);
        }
            
            
            
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        return a;
    }
    public String writeFasta(LinkedHashMap<String, DNASequence> total) {
        try(BufferedWriter bw=new BufferedWriter(new FileWriter(url+"anticuerposGrupos.fasta"));)
        {
            Iterator<String> it = total.keySet().iterator();
            while(it.hasNext())
            {
                String id=it.next();
                DNASequence idSeq=total.get(id);
                String identificador= idSeq.getOriginalHeader();
                String seq=idSeq.getSequenceAsString();
                identificador=">"+identificador;
                //Escribimos en el fichero
                bw.write(identificador);
                bw.newLine();
                bw.write(seq);
                bw.newLine();
                //Guardamos los cambios del fichero
                bw.flush();
            }            
        }catch(IOException e){
            System.out.println("Error E/S: "+e);
        }
        return url+"anticuerposGrupos";
    }
    public void writeFitness(probabiliodades p) {
        try(BufferedWriter bw=new BufferedWriter(new FileWriter(url+"grupos.txt",true));)
        {            
            //Escribimos en el fichero
            bw.write(p.getGrupo()+"-"+p.getAtriburo()+"-"+p.getProbabilidad());
            bw.newLine();
            //Guardamos los cambios del fichero
            bw.flush();
                      
        }catch(IOException e){
            System.out.println("Error E/S: "+e);
        }
    }
    public void writeAst() {
        try(BufferedWriter bw=new BufferedWriter(new FileWriter(url+"grupos.txt",true));)
        {            
            //Escribimos en el fichero
            bw.write("****");
            bw.newLine();
            //Guardamos los cambios del fichero
            bw.flush();
                      
        }catch(IOException e){
            System.out.println("Error E/S: "+e);
        }
    }
    
    public void writeProbG(ArrayList<Double> grupos,int cantG) 
    {
        try(BufferedWriter bw=new BufferedWriter(new FileWriter(url+"probGrupo.txt"));)
        {
            Iterator<Double> it = grupos.iterator();
            String ngrupo=Integer.toString(cantG);
            bw.write(ngrupo);
            bw.newLine();
        while(it.hasNext())
        {
                double numero = it.next();
                String n=Double.toString(numero);
                //Escribimos en el fichero
                bw.write(n);
                bw.newLine();
                //Guardamos los cambios del fichero
                bw.flush();
        }            
        }catch(IOException e){
            System.out.println("Error E/S: "+e);
        }
    }
    public Double ReadProbAtrb(int grupo,int atributo,int tabla) 
    {
        double probabilidad;
        String bgrupo=Integer.toString(grupo);
        String batributo=Integer.toString(atributo);
        int existe;
        int cont=0;
        int ast=0;
        try(BufferedReader br=new BufferedReader(new FileReader( url+"grupos.txt" ));)
        {
            String linea;
            while(br.ready())
            {
                if(cont>tabla)
                    break;
                if(cont==tabla)
                {
                    linea = br.readLine();
                    //System.out.println(linea);
                    if(linea.startsWith(bgrupo))
                    {
                        existe=linea.indexOf("-"+batributo+"-");
                        if(existe!=-1)
                        {                        
                            int pos=linea.lastIndexOf("-");
                            System.out.println(tabla+" "+linea+" ->"+batributo);
                            String prob=linea.substring(pos+1, linea.length());                        
                            probabilidad=Double.parseDouble(prob);
                            //System.out.println(batributo+" "+probabilidad+" "+grupo);
                            br.close();
                            return probabilidad;
                        }
                    }                   
                }
                else
                {
                    linea=br.readLine();
                }
                if(linea.startsWith("*"))
                {
                    cont++;
                }
                 
                 
            }
        }catch(IOException e){
            System.out.println("Error E/S: "+e);
        }
        return -1.0;
    }
    public double nGrupos(int opc) 
    {
        
        try(BufferedReader br=new BufferedReader(new FileReader( url+"probGrupo.txt" ));)
        {
            if(opc==0 && br.ready())
            {
                String linea=br.readLine();
                return Double.parseDouble(linea);
            }
            else
            {                
                for (int i = 0; i < opc; i++) 
                {
                    if(br.ready())
                    {
                        br.readLine();
                    }
                }
                if(br.ready())
                {
                    String linea=br.readLine();
                    return Double.parseDouble(linea);
                }
                
            }
        }catch(IOException e){
            System.out.println("Error E/S: "+e);
        }
        return 0;
    }
}



