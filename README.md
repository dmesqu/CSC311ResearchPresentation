# CSC311ResearchPresentation
Extra Credit Research Presentation Project For CSC311

![image](https://github.com/user-attachments/assets/cf3f160e-191f-4b17-a29b-bf106a0b7daa)
![image](https://github.com/user-attachments/assets/3be9e048-1804-4c91-880c-32942d7cef04)
![image](https://github.com/user-attachments/assets/21f71734-3a21-4da9-9cf6-d2c5f9fad97c)
![image](https://github.com/user-attachments/assets/90502024-1eac-4a71-9e4a-899aa3fcd738)
![image](https://github.com/user-attachments/assets/5c2400db-d9cb-45ef-a5ca-7bdc1c55ec50)
![image](https://github.com/user-attachments/assets/a83d2344-fabf-493a-8b66-ae340d8adb47)
![image](https://github.com/user-attachments/assets/ae1b2c24-bd88-4751-b63a-99954eb9111c)
![image](https://github.com/user-attachments/assets/34b34ba3-8446-4d64-b83e-4449e7a3b013)
![image](https://github.com/user-attachments/assets/af490fe2-4221-484a-aecd-c0260026018f)
````
public static void simpleCircuit() {
        Program p = new Program(2);
        Gate xGate1 = new X(0);
        Step step1 = new Step();
        step1.addGate(xGate1);
        p.addStep(step1);

        Gate hGate2 = new Hadamard(0);
        Gate xGate2 = new X(1);

        Step step2 = new Step();
        step2.addGates(hGate2, xGate2);
        p.addStep(step2);
        SimpleQuantumExecutionEnvironment sqee = new SimpleQuantumExecutionEnvironment();
        Result res = sqee.runProgram(p);
        Qubit[] qubits = res.getQubits();
        Arrays.asList(qubits).forEach(q -> System.out.println("qubit with probability on 1 = "+q.getProbability()+", measured it gives "+ q.measure()));
    }
````
![image](https://github.com/user-attachments/assets/9fa3b108-4265-413d-9f6b-308d2fff62b4)
````
public static void bellStateCreation() {
        QuantumExecutionEnvironment simulator = new SimpleQuantumExecutionEnvironment();
        Program program = new Program(2);
        Step step1 = new Step();
        step1.addGate(new Hadamard(0));
        program.addStep(step1);
        Step step2 = new Step();
        step2.addGate(new Cnot(0,1));
        program.addStep(step2);
        Result result = simulator.runProgram(program);
        Qubit[] qubits = result.getQubits();
        Qubit q0 = qubits[0];
        Qubit q1 = qubits[1];
        int v0 = q0.measure();
        int v1 = q1.measure();
        Renderer.renderProgram(program);
        Renderer.showProbabilities(program, 1000);
        }
````
![image](https://github.com/user-attachments/assets/8caf0b66-9ed1-412f-b61c-f147cc07c1c3)
````
private static void doGrover(int dim, int solution) {
        int N = 1 << dim; // Total number of possible states (2^dim)
        // Grover's algorithm requires about π/4 * sqrt(N) iterations
        double cnt = Math.PI * Math.sqrt(N) / 4;
        // Create a quantum program with 'dim' qubits
        Program p = new Program(dim);
        // Step 1: Apply Hadamard gates to all qubits to create equal superposition
        Step s0 = new Step();
        for (int i = 0; i < dim; i++) {
            s0.addGate(new Hadamard(i));
        }
        p.addStep(s0);
        // Step 2: Create Oracle and Diffusion operators
        Oracle oracle = createOracle(dim, solution); // Marks the correct solution
        oracle.setCaption("O");
        Complex[][] dif = createDiffMatrix(dim); // Diffusion (amplification) operator
        Oracle difOracle = new Oracle(dif);
        difOracle.setCaption("D");
        // Step 3: Grover iterations
        for (int i = 1; i < cnt; i++) {
            // P: Snapshot before applying Oracle
            Step s0prob = new Step("Prob " + i);
            s0prob.addGate(new ProbabilitiesGate(0));
            // O: Apply Oracle (mark solution by flipping phase)
            Step s1 = new Step("Oracle " + i);
            s1.addGate(oracle);
            // P: Snapshot after Oracle
            Step s1prob = new Step("Prob " + i);
            s1prob.addGate(new ProbabilitiesGate(0));
            // D: Apply Diffusion (amplify solution probability)
            Step s2 = new Step("Diffusion " + i);
            s2.addGate(difOracle);
            // P: Snapshot after Diffusion
            Step s3 = new Step("Prob " + i);
            s3.addGate(new ProbabilitiesGate(0));
            // Add all steps to the program
            p.addStep(s0prob);
            p.addStep(s1);
            p.addStep(s1prob);
            p.addStep(s2);
            p.addStep(s3);
        }

        System.out.println("n = " + dim + ", steps = " + cnt);
        // Run the quantum program
        Result res = sqee.runProgram(p);
        // Print out the probability of the solution at each step
        for (int i = 1; i < cnt; i++) {
            Complex[] ip0 = res.getIntermediateProbability(3 * i); // Probabilities after each diffusion
            System.out.println("results after step " + i + ": " + ip0[solution].abssqr());
        }
        System.out.println("\n");
        // Render the circuit diagram
        Renderer.renderProgram(p);
    }
    //Creates the diffusion matrix used to amplify the marked solution.
    static Complex[][] createDiffMatrix(int dim) {
        int N = 1 << dim;
        // Start with tensor product of Hadamard matrices
        Gate g = new Hadamard(0);
        Complex[][] matrix = g.getMatrix();
        Complex[][] h2 = matrix;
        for (int i = 1; i < dim; i++) {
            h2 = sqee.tensor(h2, matrix);
        }
        // Identity matrix with phase inversion for |0⟩ state
        Complex[][] I2 = new Complex[N][N];
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                I2[i][j] = (i == j) ? Complex.ONE : Complex.ZERO;
            }
        }
        I2[0][0] = Complex.ONE.mul(-1); // Flip the |0⟩ state
        // Apply H * I * H to create the diffusion matrix
        Complex[][] inter1 = mmul(h2, I2);
        Complex[][] dif = mmul(inter1, h2);
        return dif;
    }
    //Matrix multiplication utility for Complex matrices.
    static Complex[][] mmul(Complex[][] a, Complex[][] b) {
        int arow = a.length;
        int acol = a[0].length;
        int brow = b.length;
        int bcol = b[0].length;
        if (acol != brow) throw new RuntimeException("#cols a " + acol + " != #rows b " + brow);

        Complex[][] answer = new Complex[arow][bcol];
        for (int i = 0; i < arow; i++) {
            for (int j = 0; j < bcol; j++) {
                Complex el = Complex.ZERO;
                for (int k = 0; k < acol; k++) {
                    el = el.add(a[i][k].mul(b[k][j]));
                }
                answer[i][j] = el;
            }
        }
        return answer;
    }
    //Creates the Oracle which flips the phase of the solution state.
    static Oracle createOracle(int dim, int solution) {
        int N = 1 << dim;

        Complex[][] matrix = new Complex[N][N];
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (i != j) {
                    matrix[i][j] = Complex.ZERO;
                } else {
                    matrix[i][j] = (i == solution) ? Complex.ONE.mul(-1) : Complex.ONE;
                }
            }
        }
        return new Oracle(matrix);
    }
}
````
![image](https://github.com/user-attachments/assets/2f62ea55-3ab1-4c2d-ac8b-5dec894c9a2b)
![image](https://github.com/user-attachments/assets/94f9b139-8bc4-427e-8085-b655100c42d4)
![image](https://github.com/user-attachments/assets/f99bac91-7807-4a1c-bb45-67ec6f9094bd)
![image](https://github.com/user-attachments/assets/b71e4a3a-3f43-4c22-b12a-a1f361511deb)
![image](https://github.com/user-attachments/assets/6aa97174-ba96-45a0-b902-35d63cd97869)


Link to the slides: https://docs.google.com/presentation/d/141uQk7Q6I4menUY9DUsSdUjw4a_9q3tLUaw5hK8xkEk/edit?usp=sharing  

Link to API: https://github.com/redfx-quantum/strange  

Link to examples used in the demo: https://github.com/redfx-quantum/strangefx
