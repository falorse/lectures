/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   TriElem.h
 * Author: shota
 *
 * Created on 2017/07/10, 18:47
 */

#ifndef TRIELEM_H
#define TRIELEM_H

#include <Node.h>

class TriElem{
public:
    Node nodes_[3];
    
    TriElem(Node* node1,Node* node2,Node* node3){
        nodes_[0]=node1;
        nodes_[1]=node2;
        nodes_[2]=node3;
    }
};

#endif /* TRIELEM_H */

