using System;
using System.Collections;
using System.Collections.Generic;

using Rhino;
using Rhino.Geometry;

using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;



/// <summary>
/// This class will be instantiated on demand by the Script component.
/// </summary>
public class Script_Instance : GH_ScriptInstance
{
    #region Utility functions
    /// <summary>Print a String to the [Out] Parameter of the Script component.</summary>
    /// <param name="text">String to print.</param>
    private void Print(string text) { /* Implementation hidden. */ }
    /// <summary>Print a formatted String to the [Out] Parameter of the Script component.</summary>
    /// <param name="format">String format.</param>
    /// <param name="args">Formatting parameters.</param>
    private void Print(string format, params object[] args) { /* Implementation hidden. */ }
    /// <summary>Print useful information about an object instance to the [Out] Parameter of the Script component. </summary>
    /// <param name="obj">Object instance to parse.</param>
    private void Reflect(object obj) { /* Implementation hidden. */ }
    /// <summary>Print the signatures of all the overloads of a specific method to the [Out] Parameter of the Script component. </summary>
    /// <param name="obj">Object instance to parse.</param>
    private void Reflect(object obj, string method_name) { /* Implementation hidden. */ }
    #endregion

    #region Members
    /// <summary>Gets the current Rhino document.</summary>
    private readonly RhinoDoc RhinoDocument;
    /// <summary>Gets the Grasshopper document that owns this script.</summary>
    private readonly GH_Document GrasshopperDocument;
    /// <summary>Gets the Grasshopper script component that owns this script.</summary>
    private readonly IGH_Component Component;
    /// <summary>
    /// Gets the current iteration count. The first call to RunScript() is associated with Iteration==0.
    /// Any subsequent call within the same solution will increment the Iteration count.
    /// </summary>
    private readonly int Iteration;
    #endregion

    /// <summary>
    /// This procedure contains the user code. Input parameters are provided as regular arguments,
    /// Output parameters as ref arguments. You don't have to assign output parameters,
    /// they will have a default value.
    /// </summary>
    private void RunScript(Point3d atr, double damping, double grav, double dt, ref object A, ref object B, ref object C, ref object D, ref object E)
    {

        g.Move(atr, damping, grav, dt);

        A = g.GetCircles();
        B = g.GetForces();
        C = g.GetVelocities();
        //D = g.GetTrajectories();
        E = g.GetEdges();
    }

    // <Custom additional code> 

    Graph g = new Graph(20);

    class Graph
    {
        public Graph(int nodecount)
        {
            for (int i = 0; i < nodecount; ++i)
            {
                AddNode(new Point3d(i, 0.0, 5.0), Vector3d.Zero, 0.2 + i * 0.1);

            }
            for (int i = 0; i < nodecount - 1; ++i)
            {
                AddEdge(nodes[i], nodes[i + 1]);
            }
        }

        public Node AddNode(Point3d p0, Vector3d u0, double mass)
        {
            Node n = new Node(p0, u0, mass);

            nodes.Add(n);
            return n;
        }

        public Edge AddEdge(Node n0, Node n1)
        {
            Edge e = new Edge(n0, n1);
            edges.Add(e);
            return e;
        }

        public List<Node> nodes = new List<Node>();
        public List<Edge> edges = new List<Edge>();

        public void Move(Point3d atr, double damping, double g, double dt)
        {
            foreach (Node n in nodes)
            {

                if (n.pts.Count > 200)
                {
                    n.pts.RemoveAt(0);
                }
                else
                {
                    n.pts.Add(n.p);
                }

                Vector3d dv = atr - n.p;
                double dist = dv.Length;
                dv.Unitize();
                n.f = dv * Math.Exp(-0.001 * dist * dist) * 0.1 + Vector3d.ZAxis * g;

                n.Move(dt, damping);
            }

        }

        public List<Circle> GetCircles()
        {
            List<Circle> c = new List<Circle>();

            foreach (Node n in nodes)
            {
                c.Add(new Circle(n.p, 0.1 * n.m));
            }
            return c;
        }

        public List<Vector3d> GetForces()
        {
            List<Vector3d> f = new List<Vector3d>();
            foreach (Node n in nodes)
            {
                f.Add(n.f);
            }
            return f;
        }

        public List<Vector3d> GetVelocities()
        {
            List<Vector3d> v = new List<Vector3d>();
            foreach (Node n in nodes)
            {
                v.Add(n.u);
            }
            return v;
        }

        public List<Polyline> GetTrajectories()
        {
            List<Polyline> pls = new List<Polyline>();
            foreach (Node n in nodes)
            {
                Polyline pl = new Polyline(n.pts);
                pls.Add(pl);
            }
            return pls;
        }

        public List<Line> GetEdges()
        {
            List<Line> els = new List<Line>();
            foreach (Edge e in edges)
            {
                Line el = new Line(e.n0.p, e.n1.p);
                els.Add(el);
            }
            return els;
        }
    }

    //Node n = new Node(new Point3d(0.0, 0.0, 0.0), new Vector3d(0.2, -0.2, 0.0), 1.0);


    public class Node
    {
        public Node(Point3d p0, Vector3d u0, double mass)
        {
            p = p0;
            u = u0;
            m = mass;
        }

        public Point3d p = Point3d.Origin;
        public Vector3d f = Vector3d.Zero;
        public Vector3d u = Vector3d.Zero;
        public double m = 0.0;

        public List<Point3d> pts = new List<Point3d>();

        public void Move(double dt, double damping)
        {
            u *= damping;
            u += f * dt * (1 / m);
            p += u * dt;


            if (p.Z < 0.0)
            {
                p.Z = 0.0;
                if (u.Z < 0.0)
                {
                    u.Z = -u.Z;
                }
            }
        }

    }

    public class Edge
    {
        public Edge(Node _n0, Node _n1)
        {
            n0 = _n0;
            n1 = _n1;
        }
        public Node n0;
        public Node n1;
    }
    // </Custom additional code> 
}