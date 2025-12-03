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
        



        for (int i = 0; i <= m; i++) {
            for (int j = 0; j <=n; j++) {
                double val = sc.nextDouble();
                int node = i * n + j+1;
                    addEdge(this.s, node,0, val);
            }
        }
        for (int i = 0; i < m; i++) {
            for (int j = 0; j <n; j++) {
                double val = sc.nextDouble();
                int node = i * n + j+1;
                    addEdge(node,this.t,0, val);
            }
        }

        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n - 1; j++) {
                double p = sc.nextDouble();
                int u = i * n + (j+1);       
                int v = i * n + (j+2); 
                addEdge(u, v,0, p);
                addEdge(v, u,0, p);
                }
            }
        }

        for (int i = 0; i < m-1; i++) {
            for (int j = 0; j < n; j++) {
                double p = sc.nextDouble();
                int u = i * n + (j+1);       
                int v = (i+1) * n + (j+1); 
                addEdge(u, v,0, p);
                addEdge(v, u,0, p);
            }
        }
    }

}




