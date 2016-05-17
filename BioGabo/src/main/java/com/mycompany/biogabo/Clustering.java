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
        LinkedHashMap<String, DNASequence> a=io.readFasta(ruta,1);
        Iterator<String> it = a.keySet().iterator();
            while(it.hasNext())
            {                
                if(i<1)
                {
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
        GapPenalty penalty = new SimpleGapPenalty( -20, -5);
        PairwiseSequenceAligner<DNASequence, NucleotideCompound> aligner = Alignments.getPairwiseAligner(
        new DNASequence(query, AmbiguityDNACompoundSet.getDNACompoundSet()),
        new DNASequence(target, AmbiguityDNACompoundSet.getDNACompoundSet()),
        PairwiseSequenceAlignerType.GLOBAL,
        penalty, SubstitutionMatrixHelper.getNuc4_4());
        SequencePair<DNASequence, NucleotideCompound>
        alignment = aligner.getPair();
        int identical = alignment.getNumIdenticals();
        return identical / (float) target.length()>=.97 && identical / (float) query.length()>=.97;
    }
    
    public void crearGrupos(DNASequence query, DNASequence target,boolean similar) throws CompoundNotFoundException
    {
        Clustering cl= new Clustering();
        
        if(grupo.isEmpty())
        {
            if (similar)
            {
                grupo.put(query.getOriginalHeader(),query);
                query.setOriginalHeader(query.getOriginalHeader()+".G1");
                target.setOriginalHeader(target.getOriginalHeader()+".G1");
                total.put(query.getOriginalHeader(), query);
                total.put(target.getOriginalHeader(), target);
                datos.put(1,2);
                grupos++;
            }
            else
            {
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
        {
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
