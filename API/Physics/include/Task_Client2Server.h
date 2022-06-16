//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//          PPPPP  H   H  EEEEE  N    N  GGGGG  L      EEEEE  III         +
//          P   P  H   H  E      NN   N  G      L      E       I          +
//          PPPPP  HHHHH  EEEEE  N N  N  G  GG  L      EEEEE   I          +
//          P      H   H  E      N  N N  G   G  L      E       I          +
//          P      H   H  EEEEE  N    N  GGGGG  LLLLL  EEEEE  III         +
//------------------------------------------------------------------------+
//          Platform for Hybrid Engineering Simulation of Flows           +
//          China Aerodynamics Research and Development Center            +
//                     (C) Copyright, Since 2010                          +
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      Task_Client2Server.h
//! @brief     Communication pattern of send data from client to server.
//!            It is not recommended for the learners in the first stage.
//! @author    Bell, Zhang Jian, He Xin.

#pragma once
#include "Task.h"
namespace PHSPACE
{
//! @brief Class of Task_Client2Server, send data from client to server.
class LIB_EXPORT Task_Client2Server : public Task
{
public:
    Task_Client2Server();
    ~Task_Client2Server();

private:
    void PreTask(ActionKey *actkey);
    void PostTask(ActionKey *actkey);
    void MainTask(ActionKey *actkey);
    void Client2Server(ActionKey *actkey, int iZone);
};

}