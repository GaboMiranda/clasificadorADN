/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.mycompany.biogabo;

/**
 *
 * @author morelos
 */
public class probabiliodades {
    private int grupo;
    private int atriburo;
    private double probabilidad;

    public probabiliodades(int grupo, int atriburo, double probabilidad) {
        this.grupo = grupo;
        this.atriburo = atriburo;
        this.probabilidad = probabilidad;
    }

    public probabiliodades() {
        grupo = 0;
        atriburo = 0;
        probabilidad = 0;
    }

    public int getGrupo() {
        return grupo;
    }

    public int getAtriburo() {
        return atriburo;
    }

    public double getProbabilidad() {
        return probabilidad;
    }

    public void setGrupo(int grupo) {
        this.grupo = grupo;
    }

    public void setAtriburo(int atriburo) {
        this.atriburo = atriburo;
    }

    public void setProbabilidad(double probabilidad) {
        this.probabilidad = probabilidad;
    }
    
    public void mostrar(probabiliodades p){
        System.out.println("grupo: "+p.getGrupo()+"\tatributo: "+p.getAtriburo()+"\tproba: "+p.getProbabilidad());
        
    }
    
}
