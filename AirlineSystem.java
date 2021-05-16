/*
Theresa Hutchins
CS1501
Assignment 4
Apr 19, 2021
*/

import java.util.*;
import java.io.*;

/*
 Some code adapted from previous labs and Segewick's Algorithms textbook.
This program represents an Airline Database system. The airline is
represented by a graph. The graph was made by using an adjacency list containing
linked edges in a linked list. These edges contain values to and from, as well
as weights (which is the distance of the flight) as well as the cost of the trip.
*/

@SuppressWarnings("unchecked")
public class AirlineSystem {
  private String [] cityNames = null; //list of the city names 
  private Graph G = null; //graph containing all trip information
  private static Scanner scan = null;
  private static final int INFINITY = Integer.MAX_VALUE;

  /**
  * Main
  */
  public static void main(String[] args) throws IOException {
    AirlineSystem airline = new AirlineSystem(); // creates new airline system
    scan = new Scanner(System.in);
    while(true){
      switch(airline.menu()){
        case 1: //read in data from a file
          airline.readGraph();
          break;
        case 2: //print entire graph
          airline.printGraph();
          break;
        case 3: //print minimum spanning tree
          airline.minTree();
          break;
        case 4://shortest distance based on miles
          airline.shortestDistance();
          break;
        case 5: //shortest distance based on price
          airline.shortestHops();
          break;
        case 6: //shortest distacne based on hops
          airline.shortestHops();
          break;
        case 7: //print trips cheaper than user input price
          airline.lowestPrice();
          break;
        case 8: //add a route 
          airline.addRoute();
          break;
        case 9: // remove a route
          airline.removeRoute();
          break;
        case 10: //quit and save to file
          airline.exportGraph(); //writes output to file
          scan.close();
          System.exit(0);
          break;
        default:
          System.out.println("Incorrect option.");
      }
    }
  }
  /*
  * Menu 
  */
  private int menu(){
    System.out.println("********************************************");
    System.out.println("Welcome to SuperFLy Airlines!");
    System.out.println("********************************************");
    System.out.println("_█████"+
                      "\n__█__███"+
                      "\n__██___███"+
                      "\n___██____██"+
                      "\n____██_____██"+
                      "\n______██____████████████████████████████"+
                      "\n_______██_______________________________██"+
                      "\n_________██______________________________███"+
                      "\n_________██________________________██████"+
                      "\n________██____________________██████"+
                      "\n________██_______________██████"+
                      "\n________██_____________████"+
                      "\n________██______________██"+
                      "\n________██_______██_______███"+
                      "\n_______██________███________██"+
                      "\n_______██_______██__██_______███"+
                      "\n_______██______██_____██______██"+
                      "\n_______██_____██_______███_____████████████"+
                      "\n_______██_____██_________███_____█________█████"+
                      "\n______██_____██___________███____________████"+
                      "\n______██____██_____________██_____███████"+
                      "\n______██____█______________██____██"+
                      "\n______██___██_____________██____██"+
                      "\n______█___██______________██___██"+
                      "\n______██_██_______________██__██"+
                      "\n_______███_________________██_██"+
                      "\n____________________________███"+
                      "\n_____________________________█");
    System.out.println("\n********************************************");
    System.out.println("\n1. Read data from a file.");
    System.out.println("2. Display all routes.");
    System.out.println("3. Compute and display the minimum spanning tree based on distance.");
    System.out.println("4. Compute Shortest path based on number of miles.");
    System.out.println("5. Compute Shortest path based on price.");
    System.out.println("6. Compute Shortest path based on number of hops.");
    System.out.println("7. Compute Trips Cheaper than a given price.");
    System.out.println("8. Add a route to the schedule.");
      System.out.println("9. Remove a route from the schedule.");
    System.out.println("10. Exit.");
    System.out.println("********************************************");
    System.out.print("Please choose a menu option (1-10): ");
    int choice = Integer.parseInt(scan.nextLine());
    return choice;
  }
  /*
  * This method takes in a file typed in by the user.
  * It is then added edge by edge to the graph.
  */
  private void readGraph() throws IOException {
    System.out.println("Please enter graph filename:");
    String fileName = scan.nextLine();
    Scanner fileScan = new Scanner(new FileInputStream(fileName));
    int v = Integer.parseInt(fileScan.nextLine());
    G = new Graph(v);

    cityNames = new String[v]; //initialize city list to the size of total vertices
    for(int i=0; i<v; i++){
      cityNames[i] = fileScan.nextLine(); //adds each city to the list
    }

    while(fileScan.hasNextLine()){ //adds all attributes to an edge which is then added to graph G
      int from = fileScan.nextInt();
      int to = fileScan.nextInt();
      int weight = fileScan.nextInt();
      double cost = fileScan.nextDouble();
      G.addEdge(new Edge(from-1, to-1, weight, cost));
      if (!fileScan.hasNextLine()){
        break;
      }
      fileScan.nextLine();
    }
    fileScan.close();
    System.out.println("Data imported successfully.");
    System.out.print("Please press ENTER to continue ...");
    scan.nextLine();
  }
  /*
  * Method that starts finding the minimum spanning tree.
  */
  private void minTree(){ //I dont think this works
    G.kruskalMST();
  }
  /*
  * Method that writes resulting graph to an output file
  * in the same format as the input.
  */
  private void exportGraph() throws IOException {
    System.out.println("Please enter the destination file name (ending in .txt): ");
    String fileName = scan.nextLine();
    FileWriter fileWriter = new FileWriter(fileName);
    fileWriter.write(G.V + "\n");
    for(int i=0;i<G.V;i++){
      fileWriter.write(cityNames[i]);
      if(i<G.V-1) fileWriter.write("\n");
    }

    for(int i=0;i<G.V;i++){
      for (Edge e : G.adj(i)){
        fileWriter.write( "\n" + Integer.toString(e.from+1) + " ");
        fileWriter.write(Integer.toString(e.to+1) + " ");
        fileWriter.write(Integer.toString(e.weight) + " ");
        fileWriter.write(Double.toString(e.cost));
      }
    } fileWriter.close();
  }

  /*
  * This method prints out each vertex as well as all destinations vistable from
  * that vertex. It also lists the cost and distance along with each destination.
  */
  private void printGraph() {
    if(G == null){
      System.out.println("Please import a graph first (option 1).");
      System.out.print("Please press ENTER to continue ...");
      scan.nextLine();
    } else {
      for (int i = 0; i < G.V; i++) {
        System.out.print(cityNames[i] + ": ");
        for (Edge e : G.adj(i)) { //try to print when no trips are available?
          if(!e.good) System.out.print(" ");
          else {System.out.print("\n\t[Destination: " + cityNames[e.to()] + ", " + "Distance: " + e.weight + 
                                 "mi, Price: $" + e.cost + "]");
        }
       }
        System.out.println();
      }
      System.out.print("Please press ENTER to continue ...");
      scan.nextLine();

    }
  }

  /*
  * This method calculates the shortest distance between two locations using
  * dijkstras algorithm. It is based on total milage.
  */
  private void shortestDistance() {
      if(G == null){
        System.out.println("Please import a graph first (option 1).");
        System.out.print("Please press ENTER to continue ...");
        scan.nextLine();
      } else {
        for(int i=0; i<cityNames.length; i++){
          System.out.println(i+1 + ": " + cityNames[i]);
        }
        System.out.print("Please enter source city (1-" + cityNames.length + "): ");
        int source = Integer.parseInt(scan.nextLine());
        System.out.print("Please enter destination city (1-" + cityNames.length + "): ");
        int destination = Integer.parseInt(scan.nextLine());
        source--;
        destination--;
        G.dijkstras(source, destination);
        if(!G.marked[destination]){
          System.out.println("There is no route from " + cityNames[source]
                              + " to " + cityNames[destination]);
        } else {
          Stack<Integer> path = new Stack<>(); //stack of edgeTo points 
          for (int x = destination; x != source; x = G.edgeTo[x]){
              path.push(x);
          }
          System.out.print("The shortest route from " + cityNames[source] +
                             " to " + cityNames[destination] + " has " +
                             G.distTo[destination] + " miles: ");

          int prevVertex = source;
          System.out.print(cityNames[source] + " ");
          while(!path.empty()){
            int v = path.pop();
            System.out.print(G.distTo[v] - G.distTo[prevVertex] + " "
                             + cityNames[v] + " ");
            prevVertex = v;
          }
          System.out.println();

        }
        System.out.print("Please press ENTER to continue ...");
        scan.nextLine();
      }
  }

  /* 
  * addRoute() allows a user to add a route from a given vertex to a destination
  * vertex. They are also able to set the distance and price of the flight. 
  */
  private void addRoute(){
    if(G == null){
      System.out.println("Please import a graph first (option 1).");
      System.out.print("Please press ENTER to continue ...");
      scan.nextLine();
      } else {
    for(int i=0; i<cityNames.length; i++){
          System.out.println(i+1 + ": " + cityNames[i]);
        }
    System.out.print("Please enter source city (1-" + cityNames.length + "): ");
    int source = Integer.parseInt(scan.nextLine());
    System.out.print("Please enter destination city (1-" + cityNames.length + "): ");
    int destination = Integer.parseInt(scan.nextLine());
    System.out.print("Please enter the flight's distance (miles): ");
    int distance = Integer.parseInt(scan.nextLine());
    System.out.print("Please enter the flight's price (USD): ");
    int price = Integer.parseInt(scan.nextLine());
    Edge e = new Edge(source-1, destination-1, distance, price);
    G.addEdge(e);
    //System.out.println(e);
    System.out.println("Route added: from " + cityNames[source-1] + " to " + cityNames[destination-1] +
                        " with a distance of " + distance + " miles, and a price of " + price + "USD.");
     System.out.print("Please press ENTER to continue ...");
      scan.nextLine();
    }

  }

  /* 
  * removeRoute() allows a user to remove a route from a given vertex to a destination
  * vertex. 
  */
  private void removeRoute(){
    if(G == null){
      System.out.println("Please import a graph first (option 1).");
      System.out.print("Please press ENTER to continue ...");
      scan.nextLine();
      } else {
        for(int i=0; i<cityNames.length; i++){
          System.out.println(i+1 + ": " + cityNames[i]);
         }
      System.out.print("Please enter source city (1-" + cityNames.length + "): ");
      int source = Integer.parseInt(scan.nextLine());
      System.out.print("Please enter destination city (1-" + cityNames.length + "): ");
      int destination = Integer.parseInt(scan.nextLine());

      G.removeEdge(source-1, destination-1);

      System.out.println("Route removed: from " + cityNames[source-1] + " to " + cityNames[destination-1]);
      System.out.print("Please press ENTER to continue ...");
      scan.nextLine();
    }

  }

  /*
  * This method calculates the shortes route available from a source to 
  * a destination based on hops. It used the breadth first search algorithm.
  */
  private void shortestHops() {
    if(G == null){
      System.out.println("Please import a graph first (option 1).");
      System.out.print("Please press ENTER to continue ...");
      scan.nextLine();
      } else{
      for(int i=0; i<cityNames.length; i++){
        System.out.println(i+1 + ": " + cityNames[i]);
      }
      System.out.print("Please enter source city (1-" + cityNames.length + "): ");
      int source = Integer.parseInt(scan.nextLine());
      System.out.print("Please enter destination city (1-" + cityNames.length + "): ");
      int destination = Integer.parseInt(scan.nextLine());
      source--; //indexes are different so decrement
      destination--;
      G.bfs(source); //runs the graph on the bfs method
      if(!G.marked[destination]){
        System.out.println("There is no route from " + cityNames[source]
                            + " to " + cityNames[destination]);
      } else { //Show distance required and the path from the source to the destination.
        Stack<Integer> path = new Stack<>();   //Use a stack to construct the shortest path from the 
          for (int x = destination; x != source; x = G.edgeTo[x]){ // edgeTo array then print the number of hops (from the distTo array) and the path.
              path.push(x);
          }
          System.out.print("The shortest route from " + cityNames[source] +
                             " to " + cityNames[destination] + " has " +
                             G.distTo[destination] + " hop(s): ");
           int prevVertex = source;
          System.out.print(cityNames[source] + " ");
          while(!path.empty()){
            System.out.print("-> ");
            int v = path.pop();
            System.out.print(cityNames[v] + " ");
            prevVertex = v;
          }
          System.out.println();
      }
      System.out.print("Please press ENTER to continue ...");
      scan.nextLine();
    }
  }
  /*
  * Method to provide user with all flights obtainable that are beneath a given price range.
  * Right now it doesn't do connected flights, just baseline flights based on individual
  * edge costs.
  */
  private void lowestPrice() {
    if(G == null){
      System.out.println("Please import a graph first (option 1).");
      System.out.print("Please press ENTER to continue ...");
      scan.nextLine();
      } else{

      System.out.print("Please enter maximum cost (USD): ");
      int ucost = Integer.parseInt(scan.nextLine());
      System.out.print("These flights are less than $"+ ucost+ ": \n");

      for(int i=0;i<G.V;i++){
        for (Edge e : G.adj(i)){
          if(e.cost < ucost){
            System.out.println("[From "+ cityNames[e.from] + " to " + cityNames[e.to] + " costs " + e.cost + "$]");
          }
        }
      }
     
      System.out.print("Please press ENTER to continue ...");
      scan.nextLine();
    }
  }

/*
* class creating undirected graph using an adjacency list
* Some code adapted from ALgorithms by Segewick and previous labs
*/
public class Graph {
  private final int V; //vertices
  private int e;       //edges
  private LinkedList<Edge>[] adj; //adjacency list
  private boolean[] marked;  // marked[v] = is there an s-v path
  private int[] edgeTo;      // edgeTo[v] = previous edge on shortest s-v path
  private int[] distTo;      // distTo[v] = number of edges shortest s-v path
  private Edge edges[];
  private Queue<Edge> mst;

  

  public Graph(int V){
    this.V = V;
    @SuppressWarnings("unchecked")
    LinkedList<Edge>[] temp = (LinkedList<Edge>[]) new LinkedList[V];
    adj=temp;
    for(int v = 0; v < V; v++){
      adj[v] = new LinkedList<Edge>(); //initialize empty lists
    }
  }
  //returns number of vertices 
 public int V() { return V; }

  //returns number of edges
  public int E() { return e; }

  /*
  * union find class to be used with Kruskal's
  * Adapted from Segewick
  */
  public class UF { 
 
    private int[] id; //id array
    private int count;

    public UF(int N){ //creates new union find
      count = N;
      id = new int[N];
      for(int i = 0; i<N; i++) 
        id[i]=i;
    }
      public int count() { return count; } //returns count

      public boolean connected(int p, int q){ //determines if p and q are the same
        return find(p) == find(q);
      }

      public int find( int p){ //finds value of id[p]
        return id[p];
      }

      public void union( int p, int q){ //unionizes p and q
        int pID = find(p);
        int qID = find(q);

        if(pID == qID) return;

        for(int i =0; i< id.length; i++){
          if(id[i]==pID) id[i]= qID;
        count--;
        }
    }
    public int either(){ return V; } //returns V

  }
 
  /* 
  * Method to create the minimum spanning tree withing the graph.
  * I used kruscals greedy algorithm to complete this task.
  * I am still confused about this implementation so it doesn't work
  * to par.
  */
  public void kruskalMST(){ //I can't tell if I did this correctly becuase nothing prints.
    Edge[] mintree = new Edge[V-1];
    Graph N = new Graph(V-1);


    Queue<Edge> mst = new LinkedList<Edge>();
    MinPQ<Edge> pq = new MinPQ<Edge>();
    for(Edge e : adj(V-1)){
      pq.insert(e);
      UF uf = new UF(G.V);

      while(!pq.isEmpty() && mst.size() < G.V-1){
        Edge f = pq.delMin(); //get min weight edge
        int v = f.to(), w = f.other(G.V);
        if(uf.connected(v,w)) continue;

        uf.union(v,w);
        mst.add(f);
        N.addEdge(f); //adds the edge f to the new graph to br printed out
        System.out.print(cityNames[f.to]);

      }
      //try to print the resulting tree made of min edges
        for (int i = 0; i < V-1; i++) {
        System.out.print(cityNames[i] + ": ");
        for (Edge o : N.adj(i)) { 
          if(!o.good) System.out.print(" ");
          else {System.out.print("\n\t[Destination: " + cityNames[o.to()] + ", " + "Distance: " + e.weight + 
                          "mi, Price: $" + o.cost + "]");
        }
       }
        System.out.println();
      }
    }
  
    System.out.print("Please press ENTER to continue ...");
    scan.nextLine();

    
    }

  public Iterable<Edge> edges() { return mst; }
  
  //adds an edge - does work
  public void addEdge(Edge edge) {
    int from = edge.from();
    adj[from].add(edge);
    e++;
    }
  //removes an edge - doesn't work (works for showing the total routes only.)
  public void removeEdge(int from, int to) {
    int i = 0;
    Edge remove = new Edge(from, to, 0 , 0.0);
    for (i = 0; i < G.V-1; i++) {
        for (Edge e : G.adj(i)) {
          if(e.to == to && e.from == from) e.good = false;
            remove = e;

          
        }
      } 
     adj[i].remove(remove);
  }
   
  //iterable adj[]
  public Iterable<Edge> adj(int v) {
      return adj[v];
    }

  //breadth first search- does work
  //adapted from lab 10
  public void bfs(int source) {
    marked = new boolean[this.V];
    distTo = new int[this.V];
    edgeTo = new int[this.V];

    Queue<Integer> q = new LinkedList<Integer>();
    for (int i = 0; i < V; i++){
      distTo[i] = INFINITY;
      marked[i] = false;
    }
    distTo[source] = 0;
    marked[source] = true;
    q.add(source);

    while (!q.isEmpty()) {
      int v = q.remove();
      for (Edge w : adj(v)) {
        if (!marked[w.to()]) {
          edgeTo[w.to()] = v;
          distTo[w.to()] = distTo[v] + 1;
          marked[w.to()] = true;
          q.add(w.to());
          }
        }
      }
    }

  // dijkstras algorithm. This doesnt work yet.
  //adapted from lab 10.
  public void dijkstras(int source, int destination) {
    marked = new boolean[this.V];
    distTo = new int[this.V];
    edgeTo = new int[this.V];
    for (int i = 0; i < V; i++){
      distTo[i] = INFINITY;
      marked[i] = false;
    }
    distTo[source] = 0;
    marked[source] = true;
    int nMarked = 1;
    int current = source;
    while (nMarked < this.V) {
      for (Edge w : adj(current)) {
        if (distTo[current]+w.weight < distTo[w.to()]) { //if curr distnce is less than prev
          edgeTo[w.to()] = current;
          distTo[w.to()] = distTo[current] + edgeTo[w.to()];
          if(!marked[current]){
            marked[current]=true;
            nMarked++;
           }             //TODO:update edgeTo and distTo        
          }
        }
        //Find the vertex with minimim path distance
        //This can be done more effiently using a priority queue!
        int min = INFINITY;
        current=-1;

        for(int i=0; i<distTo.length; i++){
          if(marked[i])
            continue;
          if(distTo[i] < min){
            min = distTo[i];
            current = i;            
          }
        //   if(!marked[w.to()]){
        //     marked[w.to()]=true;
        //     nMarked++;
        //   }
        }//TODO: Update marked[] and nMarked. Check for disconnected graph.
      }
    }
  }

/*
* Edge class. Represents routes to and from cities. Contains the flights
* origin, destination, distance and cost. I added a boolean to determine if the
* route was deleted or not.
*/
public class Edge implements Comparable<Edge>{
  private int to;
  private int from;
  private double cost;
  private int weight;
  private boolean good;

  public Edge(int f, int t, int w, double c){ //creates an edge
    this.to = t;
    this.from = f;
    this.weight = w;
    this.cost=c;
    this.good=true;
  }

  public Edge(){}

  //comparable method based on weight (distance)
  public int compareTo(Edge c){
    return this.weight = c.weight;
  }

  public int other(int vert){
    if(vert == to) return from;
    return to;
  }

  public double getCost(){ //returns the cost
    return this.cost;
  }
  public int to(){
    return this.to;
  }
  public int from(){
    return this.from;
  }

  public void setNext(int n){
    this.to = n;
  }

  public void setPrev(int p){
    this.from = p;
  }

}

  /*
  * Below is just javas MinPQ code written by Segewick. I used the authors code 
  * for Kruskals and kept getting errors so this fixed it. WIll be inquiring with
  * Dr about this during office hours (becuase I am confused). I am certain there
  * was a better way to do this but I wanted to focus on other things.
  */
  public class MinPQ<Key> implements Iterable<Key> {
    private Key[] pq;                    // store items at indices 1 to n
    private int n;                       // number of items on priority queue
    private Comparator<Key> comparator;  // optional comparator

    /**
     * Initializes an empty priority queue with the given initial capacity.
     *
     * @param  initCapacity the initial capacity of this priority queue
     */
    public MinPQ(int initCapacity) {
        pq = (Key[]) new Object[initCapacity + 1];
        n = 0;
    }

    /**
     * Initializes an empty priority queue.
     */
    public MinPQ() {
        this(1);
    }

    /**
     * Initializes an empty priority queue with the given initial capacity,
     * using the given comparator.
     *
     * @param  initCapacity the initial capacity of this priority queue
     * @param  comparator the order in which to compare the keys
     */
    public MinPQ(int initCapacity, Comparator<Key> comparator) {
        this.comparator = comparator;
        pq = (Key[]) new Object[initCapacity + 1];
        n = 0;
    }

    /**
     * Initializes an empty priority queue using the given comparator.
     *
     * @param  comparator the order in which to compare the keys
     */
    public MinPQ(Comparator<Key> comparator) {
        this(1, comparator);
    }

    /**
     * Initializes a priority queue from the array of keys.
     * <p>
     * Takes time proportional to the number of keys, using sink-based heap construction.
     *
     * @param  keys the array of keys
     */
    public MinPQ(Key[] keys) {
        n = keys.length;
        pq = (Key[]) new Object[keys.length + 1];
        for (int i = 0; i < n; i++)
            pq[i+1] = keys[i];
        for (int k = n/2; k >= 1; k--)
            sink(k);
        assert isMinHeap();
    }

    /**
     * Returns true if this priority queue is empty.
     *
     * @return {@code true} if this priority queue is empty;
     *         {@code false} otherwise
     */
    public boolean isEmpty() {
        return n == 0;
    }

    /**
     * Returns the number of keys on this priority queue.
     *
     * @return the number of keys on this priority queue
     */
    public int size() {
        return n;
    }

    /**
     * Returns a smallest key on this priority queue.
     *
     * @return a smallest key on this priority queue
     * @throws NoSuchElementException if this priority queue is empty
     */
    public Key min() {
        if (isEmpty()) throw new NoSuchElementException("Priority queue underflow");
        return pq[1];
    }

    // resize the underlying array to have the given capacity
    private void resize(int capacity) {
        assert capacity > n;
        Key[] temp = (Key[]) new Object[capacity];
        for (int i = 1; i <= n; i++) {
            temp[i] = pq[i];
        }
        pq = temp;
    }

    /**
     * Adds a new key to this priority queue.
     *
     * @param  x the key to add to this priority queue
     */
    public void insert(Key x) {
        // double size of array if necessary
        if (n == pq.length - 1) resize(2 * pq.length);

        // add x, and percolate it up to maintain heap invariant
        pq[++n] = x;
        swim(n);
        assert isMinHeap();
    }

    /**
     * Removes and returns a smallest key on this priority queue.
     *
     * @return a smallest key on this priority queue
     * @throws NoSuchElementException if this priority queue is empty
     */
    public Key delMin() {
        if (isEmpty()) throw new NoSuchElementException("Priority queue underflow");
        Key min = pq[1];
        exch(1, n--);
        sink(1);
        pq[n+1] = null;     // to avoid loitering and help with garbage collection
        if ((n > 0) && (n == (pq.length - 1) / 4)) resize(pq.length / 2);
        assert isMinHeap();
        return min;
    }


   /***************************************************************************
    * Helper functions to restore the heap invariant.
    ***************************************************************************/

    private void swim(int k) {
        while (k > 1 && greater(k/2, k)) {
            exch(k, k/2);
            k = k/2;
        }
    }

    private void sink(int k) {
        while (2*k <= n) {
            int j = 2*k;
            if (j < n && greater(j, j+1)) j++;
            if (!greater(k, j)) break;
            exch(k, j);
            k = j;
        }
    }

   /***************************************************************************
    * Helper functions for compares and swaps.
    ***************************************************************************/
    private boolean greater(int i, int j) {
        if (comparator == null) {
            return ((Comparable<Key>) pq[i]).compareTo(pq[j]) > 0;
        }
        else {
            return comparator.compare(pq[i], pq[j]) > 0;
        }
    }

    private void exch(int i, int j) {
        Key swap = pq[i];
        pq[i] = pq[j];
        pq[j] = swap;
    }

    // is pq[1..n] a min heap?
    private boolean isMinHeap() {
        for (int i = 1; i <= n; i++) {
            if (pq[i] == null) return false;
        }
        for (int i = n+1; i < pq.length; i++) {
            if (pq[i] != null) return false;
        }
        if (pq[0] != null) return false;
        return isMinHeapOrdered(1);
    }

    // is subtree of pq[1..n] rooted at k a min heap?
    private boolean isMinHeapOrdered(int k) {
        if (k > n) return true;
        int left = 2*k;
        int right = 2*k + 1;
        if (left  <= n && greater(k, left))  return false;
        if (right <= n && greater(k, right)) return false;
        return isMinHeapOrdered(left) && isMinHeapOrdered(right);
    }


    /**
     * Returns an iterator that iterates over the keys on this priority queue
     * in ascending order.
     * <p>
     * The iterator doesn't implement {@code remove()} since it's optional.
     *
     * @return an iterator that iterates over the keys in ascending order
     */
    public Iterator<Key> iterator() {
        return new HeapIterator();
    }

    private class HeapIterator implements Iterator<Key> {
        // create a new pq
        private MinPQ<Key> copy;

        // add all items to copy of heap
        // takes linear time since already in heap order so no keys move
        public HeapIterator() {
            if (comparator == null) copy = new MinPQ<Key>(size());
            else                    copy = new MinPQ<Key>(size(), comparator);
            for (int i = 1; i <= n; i++)
                copy.insert(pq[i]);
        }

        public boolean hasNext()  { return !copy.isEmpty();                     }
        public void remove()      { throw new UnsupportedOperationException();  }

        public Key next() {
            if (!hasNext()) throw new NoSuchElementException();
            return copy.delMin();
        }
    }
    }


  

}


     




