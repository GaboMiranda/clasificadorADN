/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.mycompany.biogabo;

import java.io.IOException;
import java.util.Iterator;
import java.util.LinkedHashMap;
import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.alignment.template.GapPenalty;
import org.biojava.nbio.alignment.template.PairwiseSequenceAligner;
import org.biojava.nbio.core.alignment.matrices.SubstitutionMatrixHelper;
import org.biojava.nbio.core.alignment.template.SequencePair;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;

/**
 *
 * @author morelos
 */
public class Clustering {
    FastaIO io= new FastaIO();
    Bayesiano bayes;
    static LinkedHashMap<String, DNASequence> grupo=new LinkedHashMap<>();
    static LinkedHashMap<String, DNASequence> total=new LinkedHashMap<>();
    static int grupos=0;
    static LinkedHashMap<Integer,Integer> datos=new LinkedHashMap<>();
   
    public void cluster(String ruta,String url) throws IOException, CompoundNotFoundException
    {
        io.setUrl(url);
        Clustering cl= new Clustering();
        int i=0;
        boolean similar;
        String s11;
        String s22;
        DNASequence s1;
        DNASequence s2;
        //se realiza la llamada al metodo readFasta de la clase io enviando 
        //la ruta del archivo fasta y devuelve un objeto de tipo LinkedHashMap 
        //con un string como llave y un DNASequence que contiene todas las 
        //lecturas encontradas en el Archivo Fasta.
        LinkedHashMap<String, DNASequence> a=io.readFasta(ruta,1);
        Iterator<String> it = a.keySet().iterator();
        //se reccorre el contenido del LinkedHashMap mientras contenga 
        //un siguiente registro.
            while(it.hasNext())
            {                
                if(i<1)
                {
                    //Si es la primer iteración se realiza una segunda lectura
                    //para obtener la siguiente secuencia y poder alinear las dos secuecias
                    //si son alineadas con el porcentaje minimo requerido se agrupan 
                    //de lo contrario se cran dos grupos.
                    String query= it.next();
                    String target=it.next();
                    s1= a.get(query);
                    s2 = a.get(target);
                    s11 = s1.toString();
                    s22 = s2.toString();
                    s11=s11.substring(0, 15);
                    s22=s22.substring(0, 15);
                    similar=cl.alignment(s11, s22);
                    cl.crearGrupos(s1, s2, similar);
                    i++;
                }
                else
                {
                //En caso de que ya sea la segunda iteración ya tenemos al menos
                //un grupo para comprarar la nueva secuencia y se alinea con las
                //secuencias pertenecientes de cada uno de los grupos existentes,
                //si cumple con los parametros establecidos en el alineamiento con 
                //alguno de los grupos se queda en dicho grupo y en caso de no 
                //cumplirse con ningún grupo se crea un grupo nuevo.
                    String query = it.next();
                    s1= a.get(query);
                    s11 = s1.toString();
                    cl.crearGrupos(s1, null, true);
                }
            }
            //System.out.println(total.size());
            String archivo=io.writeFasta(total);
            bayes=new Bayesiano(archivo,grupos,total.size(),datos);
            bayes.crearTripletes();
            bayes.extraerSeq();
//            System.out.println("tam total: "+total.size()+" tam grupos: "+grupo.size());
    }
    
    public boolean alignment (String query, String target) throws CompoundNotFoundException
    {
        //método encargado de realizar un alineamiento entre dos secuencias
        GapPenalty penalty = new SimpleGapPenalty( -20, -5);
        PairwiseSequenceAligner<DNASequence, NucleotideCompound> aligner = Alignments.getPairwiseAligner(
        new DNASequence(query, AmbiguityDNACompoundSet.getDNACompoundSet()),
        new DNASequence(target, AmbiguityDNACompoundSet.getDNACompoundSet()),
        PairwiseSequenceAlignerType.GLOBAL,
        penalty, SubstitutionMatrixHelper.getNuc4_4());
        SequencePair<DNASequence, NucleotideCompound>
        alignment = aligner.getPair();
        int identical = alignment.getNumIdenticals();
        /*Devuelve un booleano que cuya condición es que la similitud entre 
        target y query sean de al menos .97 y también que la similitud entre 
                query y target sea de igual manera de minimo de .97*/
        return identical / (float) target.length()>=.97 && identical / (float) query.length()>=.97;
    }
    
    public void crearGrupos(DNASequence query, DNASequence target,boolean similar) throws CompoundNotFoundException
    {
        Clustering cl= new Clustering();
        if(grupo.isEmpty())
        { //Si grupo no está inicializado, es en caso de que se necesite crear 
            //el primer grupo o los primeros dos grupos 
            if (similar)
            {//Si las dos primeras secuncias son similares simplemente se crea un grupo
                grupo.put(query.getOriginalHeader(),query);
                query.setOriginalHeader(query.getOriginalHeader()+".G1");
                target.setOriginalHeader(target.getOriginalHeader()+".G1");
                total.put(query.getOriginalHeader(), query);
                total.put(target.getOriginalHeader(), target);
                datos.put(1,2);
                grupos++;
            }
            else
            {//En caso de que sean diferentes las dos primeras seceuncias se crean dos grupos
                grupo.put(query.getOriginalHeader(),query);
                grupo.put(target.getOriginalHeader(),target);
                query.setOriginalHeader(query.getOriginalHeader()+".G1");
                target.setOriginalHeader(target.getOriginalHeader()+".G2");
                total.put(query.getOriginalHeader(), query);
                total.put(target.getOriginalHeader(), target);
                datos.put(1,1);
                datos.put(2,1);
                grupos+=2;
            }
        }
        else
        {//Una vez que ya se ha realizado al menos el primer grupo y se alinea 
            /*con la lectura de la secuencia nueva se determina y se agrega en 
            el grupo correspondiente o se crea un nuevo gruupo en caso de no 
                    pertenecer a ningún grupo previamente creado*/
            Iterator<String> it = grupo.keySet().iterator();
            while(it.hasNext())
            {
                String id=it.next();
                DNASequence idSeq=grupo.get(id);
                String seq=idSeq.getSequenceAsString();
                seq=seq.substring(0, 15);
                similar=cl.alignment(query.getSequenceAsString().substring(0,15),seq);
                if(similar)
                {
                    String idT=idSeq.getOriginalHeader();
                    int pos=idT.indexOf(".G");
                    String cgrupo=idT.substring(pos+2, idT.length());
                    int igrupo=Integer.parseInt(cgrupo);
                    int clave=datos.get(igrupo);
                    datos.put(igrupo, clave+1);
                    query.setOriginalHeader(query.getOriginalHeader()+".G"+cgrupo);
                    total.put(query.getOriginalHeader(), query);
                    break;
                }
            }
            if(!similar)
            {
                grupo.put(query.getOriginalHeader(),query);
                query.setOriginalHeader(query.getOriginalHeader()+".G"+Integer.toString(grupos+1));                
                total.put(query.getOriginalHeader(), query);
                datos.put(grupos+1, 1);
                grupos++;
            }
        }
    }
    
    public void mostrarLista(LinkedHashMap<String, DNASequence> lista)
    {
        //método encargado de recorrer e imprimir en consola el contenido del LinkedHashMap
        Iterator<String> it = lista.keySet().iterator();
            while(it.hasNext())
            {
                String id=it.next();
                DNASequence idSeq=lista.get(id);
                String seq=idSeq.getSequenceAsString();
                System.out.println(idSeq.getOriginalHeader());
                System.out.println(seq);
            }
    }
}
