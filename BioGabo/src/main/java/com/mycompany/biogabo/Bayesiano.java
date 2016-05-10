/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.mycompany.biogabo;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedHashMap;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;

/**
 *
 * @author morelos
 */
public class Bayesiano {
    static String []combinaciones=new String[125];
    static ArrayList<ArrayList<Integer>> ListaValores = new ArrayList<>();
    static ArrayList<ArrayList<probabiliodades>> probaAtributo = new ArrayList<>();
    static ArrayList<Double> probabilidadGrupo = new ArrayList<>();
    private final String ruta;
    private final int grupos;
    private final int total;
    static ArrayList<Integer> ugrupos =  new ArrayList();
    LinkedHashMap<Integer,Integer> datos;
    FastaIO io= new FastaIO();
    Bayesiano bay;
    Bayesiano(String ruta, int grupos,int total,LinkedHashMap<Integer,Integer> datos) {
        this.ruta = ruta;
        this.grupos=grupos;
        this.total=total;
        this.datos=datos;
    }
    public void extraerSeq() throws IOException, CompoundNotFoundException
    {
        LinkedHashMap<String, DNASequence> lectura=io.readFasta(ruta,2);
        Iterator<String> it = lectura.keySet().iterator();
        while(it.hasNext())
        {
            String llave=it.next();
            DNASequence id=lectura.get(llave);
            this.compararTrip(id.getOriginalHeader(), id.getSequenceAsString());            
        }
        this.mostrarArray(ListaValores);
        this.CalcularProbabilidadGrupo();
        io.writeProbG(probabilidadGrupo,grupos);
        this.obtenerGruposUnicos();
        this.calcularProbabilidadTrip();
//        this.mostrarArrayObjetos();
//        this.mostrarArray(ListaValores);
    }
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
//        cont=125;
//        for (int i = 0; i < bases.length(); i++) 
//        {
//            for (int j = 0; j < bases.length(); j++) 
//            {
//                combinaciones[cont]=bases.charAt(i)+""+bases.charAt(j)+" ";
//                cont++;
//            }
//        }
//        cont=150;
//        for (int j = 0; j < bases.length(); j++) 
//        {
//            combinaciones[cont]=bases.charAt(j)+" "+" ";
//            cont++;
//        }
    }
    public void compararTrip(String id,String seq)
    {
        ArrayList<Integer> valores = new ArrayList<>();
        int pos=id.indexOf(".G");
        String cgrupo=id.substring(pos+2, id.length());
        int igrupo=Integer.parseInt(cgrupo);
        seq=seq.substring(0, 15);
//        int res=seq.length()%3;
//        if(res!=0)
//        {
//            res=3-res;
//            for (int i = 0; i < res; i++) 
//            {
//                seq+=" ";
//            }
//        }
        for (int i = 0; i < seq.length(); i++) 
        {
            String corte=seq.substring(i,i+3);
            for (int j = 0; j < combinaciones.length; j++) 
            {
                if(corte.compareTo(combinaciones[j])==0)
                {
//                    System.out.println(corte+" "+ combinaciones[j]+" "+j);
                    valores.add(j);
                }
            }
            i+=2;
        }
        valores.add(igrupo);
        ListaValores.add(valores);
        
    }
    public void mostrarTablaTrip()
    {
        for (int i = 0; i < combinaciones.length; i++) {
            System.out.println((i+1)+""+combinaciones[i]);
        }
    }
    
    public void mostrarArray(ArrayList<ArrayList<Integer>> Lista)
    {
        Iterator<ArrayList<Integer>> itrArrayListNumeros = Lista.iterator();
        while(itrArrayListNumeros.hasNext())
        {
                ArrayList<Integer> numeros = itrArrayListNumeros.next();
                Iterator<Integer> itrNumeros = numeros.iterator();
                while(itrNumeros.hasNext())
                {
                    int numero = itrNumeros.next();
//                    System.out.print(numero+"\t");
                }
                System.out.println();
        }
    }
    
    public void mostrarDatos()
    {
        Iterator<Integer> it = datos.keySet().iterator();
            while(it.hasNext())
            {                
                Integer key = it.next();
                System.out.println("Clave: " + key + " -> Valor: " + datos.get(key));
            }
    }
    
    public void CalcularProbabilidadGrupo()
    {
        double probabilidad;
        double a;
        double b;
        Iterator<Integer> it = datos.keySet().iterator();
            while(it.hasNext())
            {               
                Integer key = it.next();
                a = datos.get(key);
                b = total;
                probabilidadGrupo.add(a/b);
            }
    }
    
    public void calcularProbabilidadTrip(){
        int atributos=0;        
        Iterator<ArrayList<Integer>> itrArrayListNumeros = ListaValores.iterator();      
        ArrayList<Integer> numeros = itrArrayListNumeros.next();
        io.writeAst();
        for (int i = 0; i < numeros.size()-1; i++) {
            recorrerArray(i);
        }
    }
    
    public void recorrerArray(int pos)
    {
        ArrayList<Integer> atributo = new ArrayList<>();
        Iterator<ArrayList<Integer>> itrArrayListNumeros = ListaValores.iterator();
        while(itrArrayListNumeros.hasNext())
        {
            ArrayList<Integer> numeros = itrArrayListNumeros.next();
            if(!atributo.contains(numeros.get(pos)))
            {
                atributo.add(numeros.get(pos));
            }
        }
        this.comparar(atributo,pos);
    }
    
    public void comparar (ArrayList<Integer> atributo, int pos)
    {
        double suma = 0;
        double proba = 0;
        int ngrupo=0;
        int natributo=0;
        double ind;
        ArrayList<probabiliodades> valor = new ArrayList<>();
        probabiliodades p =  new probabiliodades();
        ArrayList<Integer> incidencias = new ArrayList<>();
        for (int i = 0; i < ugrupos.size(); i++) 
        {
            ngrupo = ugrupos.get(i);
            
            for (int j = 0; j < atributo.size(); j++) 
            {
                natributo = atributo.get(j);
                incidencias.add(this.incidenciasGA(natributo, ngrupo, pos));
            }
            for (int j = 0; j < incidencias.size(); j++) {
                suma+=incidencias.get(j);
            }
            for (int j = 0; j < incidencias.size(); j++) {
                p.setAtriburo(atributo.get(j));
                p.setGrupo(ngrupo);
                ind=incidencias.get(j);
                p.setProbabilidad(ind/suma);
                System.out.println(p.getGrupo()+"-"+p.getAtriburo()+"-"+p.getProbabilidad());
                io.writeFitness(p);
                //System.out.println("incidencia: "+ind+" suma: "+suma);
                //System.out.println(ngrupo+"\t"+atributo.get(j));
//                valor.add(p);
//                p.mostrar(p);
                //System.out.println("xxxxx->"+valor.get(i).getGrupo());
            }
            //System.out.println("Grupo"+valor.get(i).getGrupo());
        //probaAtributo.add(valor);
        suma=0;
        incidencias.clear();
        }
        io.writeAst();
    }       
  
    public void obtenerGruposUnicos()
    {
        int cont = 0;
        Iterator<ArrayList<Integer>> itrArrayListNumeros = ListaValores.iterator();
        while(itrArrayListNumeros.hasNext())
        {
            ArrayList<Integer> numeros = itrArrayListNumeros.next();
            if(!ugrupos.contains(numeros.get(numeros.size()-1)))
            {
                ugrupos.add(numeros.get(numeros.size()-1));
                System.out.println(numeros.get(numeros.size()-1));
                cont++;
            }
        }
    }
    
    public int incidenciasGA(int natributo,int ngrupo,int pos){
        int cont = 0;
        ArrayList<Integer> atributo = new ArrayList<>();
        Iterator<ArrayList<Integer>> itrArrayListNumeros = ListaValores.iterator();
        while(itrArrayListNumeros.hasNext())
        {
            ArrayList<Integer> numeros = itrArrayListNumeros.next();
            if(numeros.get(pos)==natributo && numeros.get(numeros.size()-1)==ngrupo)
            {
                cont++;
            }
            
        }
        return cont;
    }
    
    public void mostrarArrayObjetos()
    {
        probabiliodades p = new probabiliodades();
        Iterator<ArrayList<probabiliodades>> itrArrayListNumeros = probaAtributo.iterator();
        while(itrArrayListNumeros.hasNext())
        {
                ArrayList<probabiliodades> numeros = itrArrayListNumeros.next();
                Iterator<probabiliodades> itrNumeros = numeros.iterator();
                while(itrNumeros.hasNext())
                {
                    probabiliodades numero = itrNumeros.next();
//                    p.mostrar(numero);
                }
                System.out.println();
        }
    }
    
}