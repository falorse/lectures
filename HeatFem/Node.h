/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Node.h
 * Author: shota
 *
 * Created on 2017/07/10, 18:47
 */

#ifndef NODE_H
#define NODE_H



class Node{
public:
    VectorXY pos_;
    double t_;
    
    Node(VectorXY pos,double t){
        pos_=pos;
        t_=t;
    }
    
    
    void print(){
        
    }
};
#endif /* NODE_H */

