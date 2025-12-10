import java.util.ArrayList;
import java.util.List;

class Edge{
    int src; #source
    int des; #destination
    int dirct;# 0 sortie,1 entre
    double cap;
    double flow;
    Edge rev;


  public Edge(int src,int des,int dirct,double cap){
    this.src=src;
    this des;
    this.dirct=dirct;
    this.cap=cap;
    this.flow=0;
    this.rev=null;
  }
}
