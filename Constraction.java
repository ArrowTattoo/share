
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

class construction{
    private List<List<Edge>> reseau;

    private int m,n;
    private int s,t; //s noir,t blance
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
    public void addEdge(int src,int des,int dirct,double cap){
        Edge down = new Edge(src des,0,cap);
        Edge up = new Edge (des,src,1,0);
        down.rev=up;
        up.rev=down;
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
    

        for (int i = 0; i < m-1; i++) {
            for (int j = 0; j < n; j++) {
                double p = sc.nextDouble();
                int u = i * n + (j+1);       
                int v = (i+1) * n + (j+1); 
                addEdge(u, v,0, p);
                addEdge(v, u,0, p);
            }
        }
    
    private boolean bfs(Edge[] tmp) {
        for (int i = 0; i < numNode; i++) {
            tmp[i] = null;
        }
        //si on  peut visit c'est True,sinon False
        boolean[] visit = new boolean[numNode];
        List<Integer> l= new ArrayList<>();
        l.add(this.s);
        visit[this.s] = true;  

        while (!l.isEmpty()) {
            int u = l.remove(0); 
            for (Edge e : reseau.get(u)) {
                int v = e.des;
                if (!visit[v] && (e.cap - e.flow > 0)) {
                    tmp[v] = e;
                    visit[v] = true;
                    l.add(v);
                    if (v == t) {
                        return true;
                    }
                }
            }
        }
        return false;
    }
    public CalculFlotMAX(){

    }

    public void CalculCoupeMin(){
        //si on  peut visit c'est True,sinon False
        boolean[] visit = new boolean[numNode];
        List<Integer> l= new ArrayList<>();
        l.add(this.s);
        visit[this.s] = true;  

        while (!l.isEmpty()) {
            int u = l.remove(0); 
            for (Edge e : reseau.get(u)) {
                int v = e.des;
                if (!visit[v] && (e.cap - e.flow > 0)) {
                    tmp[v] = e;
                    visit[v] = true;
                    l.add(v);
                }
            }
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                int node = i * m + j + 1;
                if (visit[node]) {
                    System.out.print("A "); 
                } 
                else {
                    System.out.print("B ");
                }
            }
        }


    }

}



