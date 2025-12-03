import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

class construction{
    private List<List<Edge>> reseau;

    private int m,n;
    private int s,t; #s noir,t blance
    private int numNode;

    public construction(int m,int n){
        this.s=0;
        this.t=m*n+1;
        this.numNode=m*n+2
        this.reseau=new ArrayList<>(this.numNode);
        for(int i=0;i<numNode;i++){
            reseau.add(new ArrayList<>);
        }
    }
    public void constructionReseau(Scanner filename){
        

    }



}

