package comp559.lcp;

import java.util.*;

import javax.swing.JPanel;
import javax.swing.border.TitledBorder;
import javax.vecmath.Point2d;
import javax.vecmath.Vector2d;
import javax.vecmath.Vector3d;
import javax.vecmath.Vector4d;

import mintools.parameters.BooleanParameter;
import mintools.parameters.DoubleParameter;
import mintools.parameters.IntParameter;
import mintools.swing.CollapsiblePanel;
import mintools.swing.VerticalFlowPanel;
import no.uib.cipr.matrix.DenseMatrix;

/**
 * Class for detecting and resolving collisions.  Currently this class uses penalty forces between rigid bodies.
 * @author kry
 */
public class CollisionProcessor {

    private List<RigidBody> bodies;
    
    /**
     * The current contacts that resulted in the last call to process collisions
     */
    public ArrayList<Contact> contacts = new ArrayList<Contact>();
    // Create an empty hash map
    static HashMap<Tuple<Integer,Integer>, Integer> currentPairs = new HashMap<>();

    /**
     * Creates this collision processor with the provided set of bodies
     * @param bodies
     */
    public CollisionProcessor( List<RigidBody> bodies ) {
        this.bodies = bodies;
    }
    
    /** keeps track of the time used for collision detection on the last call */
    double collisionDetectTime = 0;
    
    /** keeps track of the time used to solve the LCP based velocity update on the last call */
    double collisionSolveTime = 0;

    // some variables to help with LCP solve
    private Vector3d normal_3D = new Vector3d();
    private Vector3d forward = new Vector3d(0,0,1);
    private Vector3d t1 = new Vector3d();
    private Vector3d t2 = new Vector3d();
    private Vector2d ri = new Vector2d();
    private Vector2d rj = new Vector2d();
    private DenseMatrix final_lambda;
    static DenseMatrix previous_lambda;
    private DenseMatrix J;
//    private DenseMatrix M;
    private DenseMatrix M_inverse;
    private DenseMatrix T;
    private DenseMatrix A;
    private DenseMatrix b;
    private DenseMatrix u;
    private DenseMatrix u2;
    private DenseMatrix f_ext;
    private DenseMatrix s;
    private DenseMatrix S;
    DenseMatrix lambdas;
    private int N;
    private int K;

    /**
     * Processes all collisions 
     * @param dt time step
     */
    public void processCollisions( double dt ) {
        contacts.clear();
        Contact.nextContactIndex = 0;
        
        long now = System.nanoTime();
        broadPhase(); // search for collisions
        collisionDetectTime = ( System.nanoTime() - now ) * 1e-9;
                
        if ( contacts.size() > 0  && doLCP.getValue() ) {
            now = System.nanoTime();

            N = bodies.size();
            K = contacts.size();

            double bounce = restitution.getValue();
            double mu = friction.getValue();
            // TODO: Compute velocity update with iterative solve of contact constraint matrix.

            J = new DenseMatrix(3 * K, 6 * N);
//            M = new DenseMatrix(6 * N, 6 * N);
            M_inverse = new DenseMatrix(6 * N, 6 * N);
            T = new DenseMatrix(6 * N, 3 * K);
            A = new DenseMatrix(3 * K, 3 * K);
            b = new DenseMatrix(3 * K, 1);
            u = new DenseMatrix(6 * N, 1);
            u2 = new DenseMatrix(6 * N, 1);
            f_ext = new DenseMatrix(1, 6 * N);
            s = new DenseMatrix(7 * N, 1);
            S = new DenseMatrix(7 * N, 6 * N);
            lambdas = new DenseMatrix(3 * K, iterations.getValue() + 1);
            final_lambda = new DenseMatrix(3 * K, 1);
            if (warmStart) {
                int current_k = 0;
                for (Contact contact : contacts) {
                    int contact_i = bodies.indexOf(contact.body1);
                    int contact_j = bodies.indexOf(contact.body2);

                    if (currentPairs.containsKey(new Tuple(contact_i,contact_j))) {
                        int k = currentPairs.get(new Tuple(contact_i,contact_j));
                        lambdas.set(current_k, 0, previous_lambda.get(k, 0));
                    }
                    current_k++;
                }
            }

            for (int i = 0; i < (6 * N); i += 6) {
                RigidBody body = bodies.get(i / 6);
                for (int m = i; m < (i + 3); m++) {
//                    M.set(m, m, body.massLinear);
                    M_inverse.set(m, m, body.minv);
//                    M.set(m + 3, m + 3, body.massAngular);
                    M_inverse.set(m + 3, m + 3, body.jinv);
                }
                f_ext.set(0, i, body.force.x);
                f_ext.set(0, i+1, body.force.y);
                f_ext.set(0, i+2, 0);
                u.set(i, 0, body.v.x);
                u.set(i+1, 0, body.v.y);
                u.set(i+2, 0, 0);

                // position matrix
                setUpPositionMatrix(body, (i / 6));

                Vector3d torque = new Vector3d(0, 0, body.torque);
                Vector3d omega = new Vector3d(0, 0, body.omega);
                Vector3d tmp_I = new Vector3d(0, 0, body.massAngular * body.omega);
                Vector3d tmp2_I = new Vector3d();
                tmp2_I.cross(omega, tmp_I);
                torque.sub(tmp2_I);

                f_ext.set(0, i+3, torque.x);
                f_ext.set(0, i+4, torque.y);
                f_ext.set(0, i+5, torque.z);
                u.set(i+3,0, 0);
                u.set(i+4,0, 0);
                u.set(i+5,0, body.omega);
            }

            // algorithm stuff
            getJacobian(contacts);
            getA();
            get_b(dt, bounce);
            GaussSeidel(mu);
            updateU(dt);

            advanceRigidBodies(dt);

            collisionSolveTime = (System.nanoTime() - now) * 1e-9;
        }
    }

    private void advanceRigidBodies(double dt) {
        DenseMatrix s2 = new DenseMatrix(7 * N, 1);
        for (int i = 0; i < (7 * N); i++) {
            double r = Math.floorMod(i, 7);
            double sum = 0;
            for (int j = 0; j < (6 * N); j++) {
                sum += S.get(i, j) * u.get(j,0);
            }
            double new_s = s.get(i, 0) + (dt * sum);
            s2.set(i, 0, new_s);

            if (r == 6) {
                RigidBody body = bodies.get(i / 7);

                body.v.x = u2.get((i / 7) * 6, 0);
                body.v.y = u2.get(((i / 7) * 6) + 1, 0);
                // update position
//                body.x.x = s2.get((i / 7) * 7, 0);
//                body.x.y = s2.get(((i / 7) * 7) + 1, 0);

                //update orientation
                body.omega = u2.get(((i / 7) * 6) + 5, 0);
//                body.theta = s2.get(((i / 7) * 7) + 3, 0);
//                body.updateTransformations();
            }
        }


    }

    private void updateU( double dt ) {
        DenseMatrix tmp_U = new DenseMatrix(6 * N, 1);

        for (int row_T = 0; row_T < 6 * N; row_T++) {
            double sum = 0;
            for (int col_T = 0; col_T < 3 * K; col_T++) {
                sum += T.get(row_T, col_T) * final_lambda.get(col_T, 0);
            }
            tmp_U.set(row_T, 0, sum);
        }

        DenseMatrix tmp2_U = new DenseMatrix(6 * N, 1);

        for (int i = 0; i < 6 * N; i++) {
            double sum = 0;
            for (int j = 0; j < 6 * N; j++) {
                sum += M_inverse.get(i, j) * f_ext.get(0, j);
            }
            tmp2_U.set(i, 0, sum);
        }

        tmp2_U.scale(dt);

        u2.add(u);
        u2.add(tmp_U);
        u2.add(tmp2_U);
    }

    private void GaussSeidel( double mu ) {
        int max_itr = iterations.getValue();
        DenseMatrix L = new DenseMatrix(3 * K, 3 * K);
        DenseMatrix D = new DenseMatrix(3 * K, 3 * K);
        DenseMatrix U = new DenseMatrix(3 * K, 3 * K);

        DenseMatrix lambdas_hi = new DenseMatrix(3 * K, 1);
        DenseMatrix lambdas_lo = new DenseMatrix(3 * K, 1);

        for (int i = 0; i < (3 * K); i++) {
            D.set(i, i, A.get(i,i));
            for (int j = i + 1; j < (3 * K); j++) {
                L.set(j, i, A.get(j,i));
            }
            for (int j = 0; j < i; j++) {
                U.set(j,i, A.get(j,i));
            }
        }

        for (int itr = 1; itr < max_itr; itr++) {
            List<Integer> randomized_indeces = new ArrayList<>();
            for (int i = 0; i < (3 * K); i++) {
                randomized_indeces.add(i);
            }
            if (randomizeInnerLoop.getValue()) {
                Collections.shuffle(randomized_indeces);
            }

            for (int i : randomized_indeces) {
//            for (int i = 0; i < (3 * K); i++) {
                int r = Math.floorMod(i, 3);

                double sum1 = 0;
                double sum2 = 0;

                for (int j = 0; j < i; j++) {
                    sum1 += (L.get(i,j) * lambdas.get(j,itr));
                }
                for (int j = i+1; j < (3 * K); j++) {
                    sum2 += (U.get(i,j) * lambdas.get(j,itr - 1));
                }

                double lambda = (-sum1 - sum2 - b.get(i, 0)) / D.get(i,i);

                if (r != 0) {
                    lambdas_hi.set(i, 0, mu * lambdas.get(i-r, itr));
                    lambdas_lo.set(i, 0, - lambdas_hi.get(i, 0));

                    double maxL = Math.max(lambdas.get(i, itr), lambdas_lo.get(i,0));
                    double minL = Math.max(maxL, lambdas_hi.get(i,0));
                    lambda = minL;
                }

                lambdas.set(i, itr, lambda);

                if (itr == max_itr - 1) {
                    final_lambda.set(i, 0, lambda);
                }
            }

        }

        previous_lambda = new DenseMatrix(final_lambda);


    }

    private void setUpPositionMatrix(RigidBody body, int p) {
        int i = p * 7;
        int j = p * 6;
        s.set(i, 0, body.x.x);
        s.set(i+1, 0, body.x.y);
        s.set(i+2, 0, 0);
        s.set(i+3, 0, body.theta);
        s.set(i+4, 0, 0);
        s.set(i+5, 0, 0);
        s.set(i+6, 0, 1);

        double q_s = 1;
        double q_x = 0;
        double q_y = 0;
        double q_z = body.theta;
        S.set(i, j, 1);
        S.set(i+1, j+1, 1);
        S.set(i+2, j+2, 1);
        S.set(i+3, j+3, -q_x);
        S.set(i+3, j+4, -q_y);
        S.set(i+3, j+5, -q_z);
        S.set(i+4, j+3, q_s);
        S.set(i+4, j+4, q_z);
        S.set(i+4, j+5, -q_y);
        S.set(i+5, j+3, -q_z);
        S.set(i+5, j+4, q_s);
        S.set(i+5, j+5, q_x);
        S.set(i+6, j+3, q_y);
        S.set(i+6, j+4, -q_x);
        S.set(i+6, j+5, q_s);
    }

    private void get_b(double dt, double bounce) {
        DenseMatrix tmp_B = new DenseMatrix(6 * N, 1);
        for (int i = 0; i < (6 * N); i++) {
            double sum = M_inverse.get(i, i) * f_ext.get(0, i);
            tmp_B.set(i, 0, sum);
        }
        tmp_B.scale(dt);
        DenseMatrix tmp2_B = new DenseMatrix(6 * N, 1);
        tmp2_B.add(tmp_B);
        tmp2_B.add(u);

        // add bounce vector
        DenseMatrix b_bounce = new DenseMatrix(3 * K, 1);

        for (int k = 0; k < (3 * K); k++) {
            double sum1 = 0;
            double sum2 = 0;
            for (int j = 0; j < (6 * N); j++) {
                sum1 += J.get(k, j) * tmp2_B.get(j, 0);
                sum2 += J.get(k, j) * u.get(j, 0);
            }
            b.set(k, 0, sum1);
            b_bounce.set(k, 0, bounce * sum2);
        }

//        b.add(b_bounce);
    }

    private void getA() {
        for (int k = 0; k < (3 * K); k++) {
            for (int m_i = 0; m_i < (6 * N); m_i++) {
                double sum = M_inverse.get(m_i, m_i) * J.get(k, m_i);
                T.set(m_i, k, sum);
            }
        }

        for (int col_T = 0; col_T < (3 * K); col_T++) {
            for (int row_J = 0; row_J < (3 * K); row_J++) {
                double sum = 0;
                for (int col_J = 0; col_J < (6 * N); col_J++) {
                    sum += J.get(row_J, col_J) * T.get(col_J, col_T);
                }
                A.set(row_J, col_T, sum);
            }
        }
    }

    private boolean warmStart = false;

    private void getJacobian(ArrayList<Contact> contacts) {
        currentPairs.clear();

        int i;
        int j;
        int k = 0;
        for (Contact contact : contacts) {
            i = bodies.indexOf(contact.body1);
            j = bodies.indexOf(contact.body2);

            // warm start lambda
            if (warmStart) {
                currentPairs.put(new Tuple(i, j), k);
            }

            ri.sub(contact.contactW, contact.body1.x);
            rj.sub(contact.contactW, contact.body2.x);

            normal_3D.set(contact.normal.x, contact.normal.y, 0);
            normal_3D.normalize();

            t1.cross(normal_3D, forward);
            t2.cross(normal_3D, forward);
            t1.scale(-1);

            t1.normalize();
            t2.normalize();

            // J_i linear
            J.set(k, (2*i), - normal_3D.x);
            J.set(k, (2*i) + 1, - normal_3D.y);
            J.set(k, (2*i) + 2,  - normal_3D.z);
            J.set(k + 1, (2*i) + 0, - t1.x);
            J.set(k + 1, (2*i) + 1, - t1.y);
            J.set(k + 1, (2*i) + 2, - t1.z);
            J.set(k + 2, (2*i) + 0, - t2.x);
            J.set(k + 2, (2*i) + 1, - t2.y);
            J.set(k + 2, (2*i) + 2, - t2.z);

            // J_j linear
            J.set(k, (2 * j), normal_3D.x);
            J.set(k, (2 * j) + 1, normal_3D.y);
            J.set(k, (2 * j) + 2,  normal_3D.z);
            J.set(k + 1, (2 * j) + 0, t1.x);
            J.set(k + 1, (2 * j) + 1, t1.y);
            J.set(k + 1, (2 * j) + 2, t1.z);
            J.set(k + 2, (2 * j) + 0, t2.x);
            J.set(k + 2, (2 * j) + 1, t2.y);
            J.set(k + 2, (2 * j) + 2, t2.z);


            Vector3d ri_x1 = new Vector3d(0, -ri.y, 0);
            Vector3d ri_x2 = new Vector3d(0, 0, -ri.x);
            Vector3d ri_x3 = new Vector3d(-ri.y, ri.x, 0);

            Vector3d rj_x1 = new Vector3d(0, -rj.y, 0);
            Vector3d rj_x2 = new Vector3d(0, 0, -rj.x);
            Vector3d rj_x3 = new Vector3d(-rj.y, rj.x, 0);

            Vector3d ri_n = new Vector3d(normal_3D.dot(ri_x1), normal_3D.dot(ri_x2), normal_3D.dot(ri_x3));
            Vector3d ri_t1 = new Vector3d(t1.dot(ri_x1), t1.dot(ri_x2), t1.dot(ri_x3));
            Vector3d ri_t2 = new Vector3d(t2.dot(ri_x1), t2.dot(ri_x2), t2.dot(ri_x3));
            Vector3d rj_n = new Vector3d(normal_3D.dot(rj_x1), normal_3D.dot(rj_x2), normal_3D.dot(rj_x3));
            Vector3d rj_t1 = new Vector3d(t1.dot(rj_x1), t1.dot(rj_x2), t1.dot(rj_x3));
            Vector3d rj_t2 = new Vector3d(t2.dot(rj_x1), t2.dot(rj_x2), t2.dot(rj_x3));


            // J_i angular
            J.set(k, (2*i) + 3, -ri_n.x);
            J.set(k, (2*i) + 4, -ri_n.y);
            J.set(k, (2*i) + 5, -ri_n.z);
            J.set(k + 1, (2*i) + 3, -ri_t1.x);
            J.set(k + 1, (2*i) + 4, -ri_t1.y);
            J.set(k + 1, (2*i) + 5, -ri_t1.z);
            J.set(k + 2, (2*i) + 3, -ri_t2.x);
            J.set(k + 2, (2*i) + 4, -ri_t2.y);
            J.set(k + 2, (2*i) + 5, -ri_t2.z);

            // J_j angular
            J.set(k, (2*i) + 3, -rj_n.x);
            J.set(k, (2*i) + 4, -rj_n.y);
            J.set(k, (2*i) + 5, -rj_n.z);
            J.set(k + 1, (2*i) + 3, -rj_t1.x);
            J.set(k + 1, (2*i) + 4, -rj_t1.y);
            J.set(k + 1, (2*i) + 5, -rj_t1.z);
            J.set(k + 2, (2*i) + 3, -rj_t2.x);
            J.set(k + 2, (2*i) + 4, -rj_t2.y);
            J.set(k + 2, (2*i) + 5, -rj_t2.z);

            k += 3;
        }
    }
    
    /**
     * Checks for collisions between bodies.  Note that you can optionaly implement some broad
     * phase test such as spatial hashing to reduce the n squared body-body tests.
     * Currently this does the naive n squared collision check.
     */
    private void broadPhase() {
        // Naive n squared body test.. might not be that bad for small number of bodies 
        visitID++;

        for ( RigidBody b1 : bodies ) {
            for ( RigidBody b2 : bodies ) { // not so inefficient given the continue on the next line
                if ( b1.index >= b2.index ) continue;
                if ( b1.pinned && b2.pinned ) continue;
                narrowPhase( b1, b2 );
            }
        }
    }
    
    /**
     * Checks for collision between boundary blocks on two rigid bodies.
     * @param body1
     * @param body2
     */
    private void narrowPhase( RigidBody body1, RigidBody body2 ) {
        if ( ! useBVTree.getValue() ) {
            for ( Block b1 : body1.blocks ) {
                for ( Block b2 : body2.blocks ) {
                    processCollision( body1, b1, body2, b2 );
                }
            }
        } else {
            // TODO: implement code to use hierarchical collision detection on body pairs
            tree1 = new BVNode(body1.blocks, body1);
            tree2 = new BVNode(body2.blocks, body2);

            BVtreeSearch(body1, body2, tree1, tree2);
        }
    }

    private BVNode tree1;
    private BVNode tree2;

    public void BVtreeSearch(RigidBody body1, RigidBody body2, BVNode A, BVNode B) {
        A.boundingDisc.updatecW();
        B.boundingDisc.updatecW();

        A.visitID = visitID;
        B.visitID = visitID;

        if (A.boundingDisc.intersects(B.boundingDisc)) {
            if (A.isLeaf() && B.isLeaf()) {
                processCollision(body1, A.leafBlock, body2, B.leafBlock);
            }
            else if (A.isLeaf()) {
                BVtreeSearch(body1, body2, A, B.child1);
                BVtreeSearch(body1, body2, A, B.child2);
            }
            else if (B.isLeaf()) {
                BVtreeSearch(body1, body2, A.child1, B);
                BVtreeSearch(body1, body2, A.child2, B);
            }
            else {
                BVtreeSearch(body1, body2, A.child1, B.child1);
                BVtreeSearch(body1, body2, A.child1, B.child2);
                BVtreeSearch(body1, body2, A.child2, B.child1);
                BVtreeSearch(body1, body2, A.child2, B.child2);
            }
        }
        else {
            return;
        }
    }
    
    /** 
     * The visitID is used to tag boundary volumes that are visited in 
     * a given time step.  Marking boundary volume nodes as visited during
     * a time step allows for a visualization of those used, but it can also
     * be used to more efficiently update the centeres of bounding volumes
     * (i.e., call a BVNode's updatecW method at most once on any given timestep)
     */
    int visitID = 0;
    
    /**
     * Resets the state of the collision processor by clearing all
     * currently identified contacts, and reseting the visitID for
     * tracking the bounding volumes used
     */
    public void reset() {
        contacts.clear();
        Contact.nextContactIndex = 0;
        visitID = 0;            
    }
    
    // some working variables for processing collisions
    private Point2d tmp1 = new Point2d();
    private Point2d tmp2 = new Point2d();
    private Point2d contactW = new Point2d();
    private Vector2d force = new Vector2d();
    private Vector2d contactV1 = new Vector2d();
    private Vector2d contactV2 = new Vector2d();
    private Vector2d relativeVelocity = new Vector2d();
    private Vector2d normal = new Vector2d();
        
    /**
     * Processes a collision between two bodies for two given blocks that are colliding.
     * Currently this implements a penalty force
     * @param body1
     * @param b1
     * @param body2
     * @param b2
     */
    private void processCollision( RigidBody body1, Block b1, RigidBody body2, Block b2 ) {        
        double k = contactSpringStiffness.getValue();
        double c1 = contactSpringDamping.getValue();
        double threshold = separationVelocityThreshold.getValue();
        boolean useSpring = enableContactSpring.getValue();
        boolean useDamping = enableContactDamping.getValue();
        
        body1.transformB2W.transform( b1.pB, tmp1 );
        body2.transformB2W.transform( b2.pB, tmp2 );
        double distance = tmp1.distance(tmp2);
        if ( distance < Block.radius * 2 ) {
            // contact point at halfway between points 
            // NOTE: this assumes that the two blocks have the same radius!
            contactW.interpolate( tmp1, tmp2, .5 );
            // contact normal
            normal.sub( tmp2, tmp1 );
            normal.normalize();
            // create the contact
            Contact contact = new Contact( body1, body2, contactW, normal);
            // simple option... add to contact list...
            contacts.add( contact );
            if ( ! doLCP.getValue() ) {
                // compute relative body velocity at contact point
                body1.getSpatialVelocity( contactW, contactV1 );
                body2.getSpatialVelocity( contactW, contactV2 );
                relativeVelocity.sub( contactV1, contactV2 );
                if ( -relativeVelocity.dot( normal ) < threshold ) {
                    if ( useSpring ) {
                        // spring force
                        double interpenetration = distance - Block.radius * 2; // a negative quantity
                        force.scale( -interpenetration * k, normal );
                        body2.applyContactForceW(contactW, force);
                        force.scale(-1);
                        body1.applyContactForceW(contactW, force);
                    }
                    if ( useDamping ) {
                        // spring damping forces!
                        // vertical
                        force.scale( relativeVelocity.dot(normal) * c1, normal );                    
                        body2.applyContactForceW( contactW, force );
                        force.scale(-1);
                        body1.applyContactForceW( contactW, force );
                    }
                }
            }
        }
    }
   
    /** Stiffness of the contact penalty spring */
    private DoubleParameter contactSpringStiffness = new DoubleParameter("penalty contact stiffness", 3e3, 1, 1e5 );
    
    /** Viscous damping coefficient for the contact penalty spring */
    private DoubleParameter contactSpringDamping = new DoubleParameter("penalty contact damping", 1, 1, 1e4 );
    
    /** Threshold for the relative velocity in the normal direction, for determining if spring force will be applied. */
    private DoubleParameter separationVelocityThreshold = new DoubleParameter( "penalty separation velocity threshold (controls bounce)", 1e-9, 1e-9, 1e3 );
    
    /** Enables the contact penalty spring */
    private BooleanParameter enableContactSpring = new BooleanParameter("enable penalty contact spring", true );
    
    /** Enables damping of the contact penalty spring */
    private BooleanParameter enableContactDamping = new BooleanParameter("enable penalty contact damping", true );
    
    /** Restitution parameter for contact constraints */
    public DoubleParameter restitution = new DoubleParameter( "restitution (bounce)", 0.5, 0, 1 );
    
    /** Coulomb friction coefficient for contact constraint */
    public DoubleParameter friction = new DoubleParameter("Coulomb friction", 0.33, 0, 2 );
    
    /** Number of iterations to use in projected Gauss Seidel solve */
    public IntParameter iterations = new IntParameter("iterations for GS solve", 10, 1, 500);
    
    /** Flag for switching between penalty based contact and contact constraints */
    private BooleanParameter randomizeInnerLoop = new BooleanParameter( "Randomize inner loop of GS", true );

    /** Flag for switching between penalty based contact and contact constraints */
    private BooleanParameter doLCP = new BooleanParameter( "do LCP solve", false );

    /** Flag for enabling the use of hierarchical collision detection for body pairs */
    private BooleanParameter useBVTree = new BooleanParameter( "use BVTree", true );
    
    /**
     * @return controls for the collision processor
     */
    public JPanel getControls() {
        VerticalFlowPanel vfp = new VerticalFlowPanel();
        vfp.setBorder( new TitledBorder("Collision Processing Controls") );
        vfp.add( useBVTree.getControls() );
        vfp.add( doLCP.getControls() );
        vfp.add( iterations.getSliderControls() );
        vfp.add( randomizeInnerLoop.getControls() );
        vfp.add( restitution.getSliderControls(false) );
        vfp.add( friction.getSliderControls(false) );
        
        VerticalFlowPanel vfp2 = new VerticalFlowPanel();
        vfp2.setBorder( new TitledBorder("penalty method controls") );
        vfp2.add( contactSpringStiffness.getSliderControls(true) );
        vfp2.add( contactSpringDamping.getSliderControls(true) );
        vfp2.add( separationVelocityThreshold.getSliderControls( true ) );
        vfp2.add( enableContactDamping.getControls() );
        vfp2.add( enableContactSpring.getControls() );
        
        CollapsiblePanel cp = new CollapsiblePanel(vfp2.getPanel());
        cp.collapse();
        vfp.add( cp );        
        return vfp.getPanel();
    }
    
}
