/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Driver.h
 * Author: shota
 *
 * Created on 2017/07/10, 19:02
 */

#ifndef DRIVER_H
#define DRIVER_H

#include <vector>


class Driver{
public:
    std::vector<Node *> nodes_;
    std::vector<TriElem *> elems_;
    
    Driver(){
    }
    
    void addNode(Node* node){
        nodes_.push_back(node);
    }
            
    void addElem(TriElem* elem){
        elems_.push_back(elem);
    }
};

#endif /* DRIVER_H */

