import streamlit as st
import subprocess
import os
import numpy as np
import plotly.graph_objects as go
import sqlite3
from pathlib import Path
import shutil
import glob
from datetime import datetime
import hashlib

# -------------------------
# Constants and metadata
# -------------------------
dict_l = {0: "s", 1: "p", 2: "d", 3: "f", 4: "g", 5: "h"}
label = {'0': 's', '1': 'p', '2': 'd', '3': 'f', '4': 'g'}
element = {
    '1': 'Hydrogen', '2': 'Helium', '3': 'Lithium', '4': 'Beryllium', '5': 'Boron',
    '6': 'Carbon', '7': 'Nitrogen', '8': 'Oxygen', '9': 'Fluorine', '10': 'Neon',
    '11': 'Sodium', '12': 'Potassium'
}
CHOICE_ORDER = ["time", "body", "potential", "perturbation", "perturbation_type", "problem"]
CHOICE_PRETTY = {
    "time": lambda v: v,
    "body": lambda v: v,
    "potential": lambda v: v,
    "perturbation": lambda v: f"Perturbation: {v}",
    "perturbation_type": lambda v: v,
    "problem": lambda v: v,
}

DB_PATH = "data/simulations.db"

# -------------------------
# Initialize DB
# -------------------------
def init_db():
    try:
        Path("data").mkdir(exist_ok=True)
        conn = sqlite3.connect(DB_PATH)
        c = conn.cursor()
        c.execute('''
            CREATE TABLE IF NOT EXISTS simulations (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                username TEXT,
                problem_name TEXT,
                output_path TEXT,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        ''')
        conn.commit()
    except sqlite3.Error as e:
        st.error(f"Database error: {e}")
    finally:
        try:
            conn.close()
        except:
            pass

init_db()

# -------------------------
# Simple authentication (replace streamlit_authenticator)
# -------------------------
# Example: produce a hashed password with:
# import hashlib
# hashlib.sha256("mypassword".encode()).hexdigest()
# Place the hex digest below for the users you want to allow.
CREDENTIALS = {
    "Bose": {
        "name": "Bose",
        # Replace the value below with the SHA-256 hex digest of the real password
        "password_hash": "7676aaafb027c825bd9abab78b234070e702752f625b752e55e55b48e607e358"
    },
    # Add more users as needed, same structure
}

def sha256_hash(text: str):
    return hashlib.sha256(text.encode("utf-8")).hexdigest()

def verify_password(username: str, password: str) -> bool:
    if username not in CREDENTIALS:
        return False
    return sha256_hash(password) == CREDENTIALS[username]["password_hash"]

# -------------------------
# Streamlit configuration & session state defaults
# -------------------------
st.set_page_config(page_title="Aliah University — Atomic Physics App", layout="wide")

if "user_choice" not in st.session_state:
    st.session_state.user_choice = {k: None for k in CHOICE_ORDER}
if "last_problem" not in st.session_state:
    st.session_state.last_problem = None
if "l" not in st.session_state:
    st.session_state.l = None
if "authenticated" not in st.session_state:
    st.session_state.authenticated = False
if "username" not in st.session_state:
    st.session_state.username = None
if "output_folder" not in st.session_state:
    st.session_state.output_folder = None
# dedicated page state (was previously handled by set_choice("page", ...))
if "page" not in st.session_state:
    st.session_state.page = "home"

# -------------------------
# Helpers
# -------------------------
def breadcrumb_text():
    parts = ["Home"]
    for k in CHOICE_ORDER:
        v = st.session_state.user_choice.get(k)
        if v:
            parts.append(CHOICE_PRETTY.get(k, lambda x: x)(v))
    if st.session_state.last_problem:
        parts.append(CHOICE_PRETTY.get("problem")(st.session_state.last_problem))
    return " → ".join(parts)

def set_choice(key, value):
    # If key is 'page', set the dedicated page state and return
    if key == "page":
        st.session_state.page = value
        return
    # Otherwise set in user_choice, clear deeper choices
    st.session_state.user_choice[key] = value
    if key in CHOICE_ORDER:
        idx = CHOICE_ORDER.index(key)
        for deeper in CHOICE_ORDER[idx + 1:]:
            st.session_state.user_choice[deeper] = None
        st.session_state.last_problem = None

def remove_last_choice():
    if st.session_state.last_problem:
        st.session_state.last_problem = None
    else:
        for k in reversed(CHOICE_ORDER):
            if st.session_state.user_choice.get(k):
                st.session_state.user_choice[k] = None
                break

# -------------------------
# Pages
# -------------------------
def login_page():
    st.title("Login")
    st.markdown("**Department of Physics — Aliah University**")
    with st.form("login_form"):
        username = st.text_input("Username")
        password = st.text_input("Password", type="password")
        submitted = st.form_submit_button("Login")
        if submitted:
            if verify_password(username, password):
                st.success("Login successful")
                st.session_state.authenticated = True
                st.session_state.username = username
                st.session_state.page = "home"
            else:
                st.error("Username/password is incorrect")
    st.markdown("*Developed by Md Sahin Ahamed, Koustav Das Chakladar, Aliah University*", unsafe_allow_html=True)

def home_page():
    st.title("Atomic Physics Simulator")
    st.markdown("""
        **Welcome, Quantum Explorer!**

        Embark on a journey through the quantum realm —  
        simulate one-body, two-body, and many-body systems,  
        explore potential wells and barriers,  
        and uncover the mysteries of atomic physics.

        **Let’s dive into the quantum world!**
    """)
    if st.button("Start Simulation"):
        set_choice("page", "time_selection")
    st.markdown("---")
    st.markdown("*Developed by Md Sahin Ahamed, Koustav Das Chakladar, Aliah University*", unsafe_allow_html=True)

def time_selection_page():
    st.subheader("Select Problem Type")
    for opt in ["Time Dependent", "Time Independent"]:
        if st.button(opt):
            set_choice("time", opt)
            set_choice("page", "body_selection")
    if st.button("Back"):
        set_choice("page", "home")
    st.markdown("---")
    st.markdown("*Developed by Md Sahin Ahamed, Koustav Das Chakladar, Aliah University*", unsafe_allow_html=True)

def body_selection_page():
    st.subheader("Select System Type")
    for opt in ["One Body", "Two Body", "Three Body", "Many Body"]:
        if st.button(opt):
            set_choice("body", opt)
            set_choice("page", "potential_selection")
    if st.button("Back"):
        remove_last_choice()
        set_choice("page", "time_selection")
    st.markdown("---")
    st.markdown("*Developed by Md Sahin Ahamed, Koustav Das Chakladar, Aliah University*", unsafe_allow_html=True)

def potential_selection_page():
    st.subheader("Select Confinement Type")
    for opt in ["Potential Well", "Potential Barrier", "penetrable spatial confinement", "None"]:
        if st.button(opt):
            set_choice("potential", opt)
            set_choice("page", "perturbation_selection")
    if st.button("Back"):
        remove_last_choice()
        set_choice("page", "body_selection")
    st.markdown("---")
    st.markdown("*Developed by Md Sahin Ahamed, Koustav Das Chakladar, Aliah University*", unsafe_allow_html=True)

def perturbation_selection_page():
    st.subheader("Is there any Perturbation?")
    if st.button("No"):
        set_choice("perturbation", "No")
        set_choice("page", "problem_selection")
    if st.button("Yes"):
        set_choice("perturbation", "Yes")
        set_choice("page", "perturbation_type_selection")
    if st.button("Back"):
        remove_last_choice()
        set_choice("page", "potential_selection")
    st.markdown("---")
    st.markdown("*Developed by Md Sahin Ahamed, Koustav Das Chakladar, Aliah University*", unsafe_allow_html=True)

def perturbation_type_selection_page():
    st.subheader("Select Perturbation Type")
    for opt in ["Electric Field Perturbation", "Magnetic Field Perturbation", "Harmonic Perturbation"]:
        if st.button(opt):
            set_choice("perturbation_type", opt)
            set_choice("page", "problem_selection")
    if st.button("Back"):
        remove_last_choice()
        set_choice("page", "perturbation_selection")
    st.markdown("---")
    st.markdown("*Developed by Md Sahin Ahamed, Koustav Das Chakladar, Aliah University*", unsafe_allow_html=True)

def problem_selection_page():
    st.subheader("Available Problems")
    problems = []
    if (st.session_state.user_choice["time"] == "Time Independent" and
        st.session_state.user_choice["body"] == "Two Body" and
        st.session_state.user_choice["potential"] == "Potential Well" and
        st.session_state.user_choice["perturbation"] == "No"):
        problems = [("Hydrogen atom in an attractive cage", "h-atom-in-att-cage/run.jl")]
    elif (st.session_state.user_choice["time"] == "Time Independent" and
          st.session_state.user_choice["body"] == "Two Body" and
          st.session_state.user_choice["potential"] == "penetrable spatial confinement" and
          st.session_state.user_choice["perturbation"] == "No"):
        problems = [("Free Two-Body system under penetrable spatial confinement", "Code_LMM/main.py")]
    else:
        problems = [("No problem available for this configuration", None)]
    for name, _ in problems:
        if st.button(name):
            st.session_state.last_problem = name
            set_choice("page", "input_page")
    back_target = "perturbation_type_selection" if st.session_state.user_choice.get("perturbation_type") else "perturbation_selection"
    if st.button("Back"):
        remove_last_choice()
        set_choice("page", back_target)
    st.markdown("---")
    st.markdown("*Developed by Md Sahin Ahamed, Koustav Das Chakladar, Aliah University*", unsafe_allow_html=True)

def input_page():
    st.subheader(st.session_state.last_problem or "No problem selected")
    inputs = {}
    if st.session_state.last_problem == "Hydrogen atom in an attractive cage":
        fields = [("n", "2"), ("l", "0"), ("N", "41"), ("maxR", "500"), ("v0 values", "-0.4,-0.2,-0.125,-0.08,-0.05")]
        for name, default in fields:
            inputs[name] = st.text_input(f"{name}:", value=default, key=f"input_{name}")
        if st.button("Run Simulation"):
            run_hydrogen_simulation(inputs)
    elif st.session_state.last_problem == "Free Two-Body system under penetrable spatial confinement":
        fields = [("n0", "1"), ("Z", "1"), ("l", "0"), ("V0", "4.0"), ("D0", "2.0"),
                  ("x0", "3.0"), ("V1", "4.0"), ("D1", "2.0"), ("x1", "6.0,7.0"), ("pk", "4")]
        for name, default in fields:
            inputs[name] = st.text_input(f"{name}:", value=default, key=f"input_{name}")
        if st.button("Run Simulation"):
            run_two_body_simulation(inputs)
    else:
        st.write("No input form for this problem.")
    col1, col2, col3 = st.columns([1, 1, 1])
    with col1:
        if st.button("Analyze Data"):
            set_choice("page", "plot_page")
    with col2:
        if st.button("Visualize Theory"):
            set_choice("page", "theory_plot_page")
    with col3:
        if st.button("Back"):
            remove_last_choice()
            set_choice("page", "problem_selection")
    st.markdown("---")
    st.markdown("*Developed by Md Sahin Ahamed, Koustav Das Chakladar, Aliah University*", unsafe_allow_html=True)

# -------------------------
# Simulation runners
# -------------------------
def run_hydrogen_simulation(inputs):
    try:
        n = int(inputs["n"])
        l = int(inputs["l"])
        st.session_state.l = l
        N = int(inputs["N"])
        maxR = int(inputs["maxR"])
        v0_values = [float(v.strip()) for v in inputs["v0 values"].split(",") if v.strip()]
        src_folder = f"{n}{dict_l.get(l, '')}"
        output_base = f"data/{src_folder}_{st.session_state.username}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        Path(output_base).mkdir(parents=True, exist_ok=True)
        for v0 in v0_values:
            st.write(f"Running simulation with v0 = {v0}...")
            output_dir = f"{output_base}/V={v0}"
            Path(output_dir).mkdir(parents=True, exist_ok=True)
            # Run Julia backend (if available)
            try:
                subprocess.run(
                    ["julia", "backend/h-atom-in-att-cage/run.jl", str(v0), str(n), str(l), str(N), str(maxR)],
                    check=True, cwd=os.getcwd()
                )
            except subprocess.CalledProcessError as se:
                st.warning(f"Julia run failed or not installed: {se}")
            # Run tabulator python script
            try:
                subprocess.run(
                    ["python3", "backend/h-atom-in-att-cage/tabulator.py", str(v0), str(n), str(l), str(N), str(maxR)],
                    check=True, cwd=os.getcwd()
                )
            except subprocess.CalledProcessError as se:
                st.warning(f"Tabulator script failed: {se}")
            # Move outputs (if any)
            for f in glob.glob(f"{src_folder}/*"):
                if os.path.exists(f):
                    try:
                        shutil.move(f, output_dir)
                    except Exception as e:
                        st.warning(f"Failed moving {f}: {e}")
            # Store metadata in DB
            try:
                conn = sqlite3.connect(DB_PATH)
                c = conn.cursor()
                c.execute("INSERT INTO simulations (username, problem_name, output_path) VALUES (?, ?, ?)",
                          (st.session_state.username, st.session_state.last_problem, os.path.abspath(output_dir)))
                conn.commit()
            finally:
                try:
                    conn.close()
                except:
                    pass
        st.session_state.output_folder = os.path.abspath(output_base)
        st.success("Simulation completed successfully!")
    except Exception as e:
        st.error(f"Error running hydrogen simulation: {str(e)}")

def run_two_body_simulation(inputs):
    try:
        n0 = int(inputs["n0"])
        Z = int(inputs["Z"])
        l = int(inputs["l"])
        st.session_state.l = l
        V0 = float(inputs["V0"])
        D0 = float(inputs["D0"])
        x0_values = [float(x0.strip()) for x0 in inputs["x0"].split(",") if x0.strip()]
        V1 = float(inputs["V1"])
        D1 = float(inputs["D1"])
        x1_values = [float(x1.strip()) for x1 in inputs["x1"].split(",") if x1.strip()]
        pk = int(inputs["pk"])
        src_folder = f"l={l}"
        output_base = f"data/{src_folder}_{st.session_state.username}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        Path(output_base).mkdir(parents=True, exist_ok=True)
        for x0 in x0_values:
            for x1 in x1_values:
                st.write(f"Running simulation with x0 = {x0}, x1 = {x1}...")
                output_dir = f"{output_base}/V0={V0}_x0={x0}_x1={x1}"
                Path(output_dir).mkdir(parents=True, exist_ok=True)
                try:
                    subprocess.run(
                        ["python3", "backend/Code_LMM/main.py", str(n0), str(Z), str(l), str(V0), str(D0), str(x0), str(V1), str(D1), str(x1), str(pk)],
                        check=True, cwd=os.getcwd()
                    )
                except subprocess.CalledProcessError as se:
                    st.warning(f"Two-body script failed: {se}")
                for f in glob.glob(f"{src_folder}/*"):
                    if os.path.exists(f):
                        try:
                            shutil.move(f, output_dir)
                        except Exception as e:
                            st.warning(f"Failed moving {f}: {e}")
                try:
                    conn = sqlite3.connect(DB_PATH)
                    c = conn.cursor()
                    c.execute("INSERT INTO simulations (username, problem_name, output_path) VALUES (?, ?, ?)",
                              (st.session_state.username, st.session_state.last_problem, os.path.abspath(output_dir)))
                    conn.commit()
                finally:
                    try:
                        conn.close()
                    except:
                        pass
        st.session_state.output_folder = os.path.abspath(output_base)
        st.success("Simulation completed successfully!")
    except Exception as e:
        st.error(f"Error running two-body simulation: {str(e)}")

# -------------------------
# Plot and theory pages
# -------------------------
def plot_page():
    st.subheader(f"Plotter — {st.session_state.last_problem}")
    if not st.session_state.output_folder or not os.path.exists(st.session_state.output_folder):
        st.warning("No output data found. Run a simulation first.")
        if st.button("Back"):
            set_choice("page", "input_page")
        return
    orbital = os.path.basename(st.session_state.output_folder).split('_')[0]
    subfolders = sorted([d for d in os.listdir(st.session_state.output_folder) if os.path.isdir(os.path.join(st.session_state.output_folder, d))])
    config = st.selectbox("Configuration:", subfolders)
    if config:
        sub_path = os.path.join(st.session_state.output_folder, config)
        if st.session_state.last_problem == "Hydrogen atom in an attractive cage":
            barrier_file = os.path.join(sub_path, f"Barrier_{config}.dat")
            wave_file = os.path.join(sub_path, f"one_particle_density_position_{config}.dat")
            fig = go.Figure()
            if os.path.exists(barrier_file):
                data = np.loadtxt(barrier_file)
                r, barrier = data[:, 0], data[:, 1]
                l = st.session_state.l or 0
                potential = np.where(r != 0, -1 / r + l * (l + 1) / (2 * r ** 2), 0)
                fig.add_trace(go.Scatter(x=r, y=potential, mode='lines', name='Effective Potential'))
                fig.add_trace(go.Scatter(x=r, y=barrier, mode='lines', name='Potential Well'))
                if os.path.exists(wave_file):
                    wave_data = np.loadtxt(wave_file)
                    wave = wave_data[:, 1] / np.max(np.abs(wave_data[:, 1])) * 0.3
                    fig.add_trace(go.Scatter(x=r, y=wave, mode='lines', name='Wavefunction', line=dict(dash='dash')))
                fig.update_layout(title="Potential and Wavefunction", xaxis_title="R", yaxis_title="V(R)", showlegend=True)
                st.plotly_chart(fig, use_container_width=True)
            # Tabulated data
            v0_str = config.split('V=')[1].split(',')[0] if 'V=' in config else config
            v0_clean = v0_str.replace('-', '')
            tab_file = os.path.join(os.path.dirname(st.session_state.output_folder), f"data/{orbital}_v0{v0_clean}.dat")
            if os.path.exists(tab_file):
                data = np.loadtxt(tab_file, skiprows=1)
                column_names = [
                    'R', 'Total Energy', 'Kinetic Energy', 'Potential Energy', 'Centrifugal Kinetic Energy',
                    'Expectation value of r', 'Expectation value of r^2', 'Expectation value of r^3'
                ]
                col1 = st.selectbox("X-Axis:", column_names, index=0, key="plot_x_axis")
                col2 = st.selectbox("Y-Axis:", column_names, index=1, key="plot_y_axis")
                xi = column_names.index(col1)
                yi = column_names.index(col2)
                fig = go.Figure()
                fig.add_trace(go.Scatter(x=data[:, xi], y=data[:, yi], mode='lines+markers', name=col2))
                fig.update_layout(title="Tabulated Data", xaxis_title=col1, yaxis_title=col2, showlegend=True)
                st.plotly_chart(fig, use_container_width=True)
        else:
            st.write("Plotting not implemented for this problem yet.")
    if st.button("Back"):
        set_choice("page", "input_page")
    st.markdown("---")
    st.markdown("*Developed by Md Sahin Ahamed, Koustav Das Chakladar, Aliah University*", unsafe_allow_html=True)

def theory_plot_page():
    st.subheader("Interactive Theory Plot")
    num_functions = st.selectbox("Number of Functions to Plot:", [1, 2, 3, 4, 5], index=0, key="num_functions")
    x_range = st.text_input("X-Axis Range (min,max):", value="0,50", key="x_range")
    y_range = st.text_input("Y-Axis Range (min,max):", value="-1,1" if st.session_state.last_problem == "Hydrogen atom in an attractive cage" else "-10,10", key="y_range")
    functions = []
    for i in range(num_functions):
        st.write(f"**Function {i+1}**")
        default_expr = "-v0 * ((x >= r) & (x <= r + d))" if st.session_state.last_problem == "Hydrogen atom in an attractive cage" else "v0 * exp(-((x - x0)**2)/(2*d0**2))"
        expr = st.text_input(f"Expression {i+1} (use x as variable):", value=default_expr, key=f"expr_{i}")
        default_params = "v0,r,d" if st.session_state.last_problem == "Hydrogen atom in an attractive cage" else "v0,d0,x0"
        params = st.text_input(f"Parameters {i+1} (comma separated):", value=default_params, key=f"params_{i}")
        default_values = "0.4,10.0,5.0" if st.session_state.last_problem == "Hydrogen atom in an attractive cage" else "4.0,2.0,3.0"
        values = st.text_input(f"Initial Values {i+1} (comma separated):", value=default_values, key=f"values_{i}")
        functions.append((expr, params, values))
    try:
        x_min, x_max = map(float, x_range.split(","))
        y_min, y_max = map(float, y_range.split(","))
        x = np.linspace(x_min, x_max, 1000)
        fig = go.Figure()
        colors = ['#00FF66', '#FF7A00', '#2979FF', '#FF5555', '#AA00FF']
        for i, (expr, param_str, val_str) in enumerate(functions):
            params = [p.strip() for p in param_str.split(",") if p.strip()]
            values = [float(v.strip()) for v in val_str.split(",") if v.strip()]
            if len(params) != len(values):
                st.error(f"Function {i+1}: Number of parameters and values must match.")
                return
            def func(x_arr, params, values):
                local_dict = {params[j]: values[j] for j in range(len(params))}
                local_dict.update({"x": x_arr, "exp": np.exp, "sin": np.sin, "cos": np.cos, "tan": np.tan, "log": np.log, "sqrt": np.sqrt, "pi": np.pi, "e": np.e, "abs": np.abs})
                return eval(expr, {"__builtins__": None}, local_dict)
            y = func(x, params, values)
            fig.add_trace(go.Scatter(x=x, y=y, mode='lines', name=f"Function {i+1}"))
        fig.update_layout(title="Interactive Function Plot", xaxis_title="x", yaxis_title="f(x)", showlegend=True, yaxis=dict(range=[y_min, y_max]))
        st.plotly_chart(fig, use_container_width=True)
    except Exception as e:
        st.error(f"Error generating plot: {str(e)}")
    if st.button("Back"):
        set_choice("page", "input_page")
    st.markdown("---")
    st.markdown("*Developed by Md Sahin Ahamed, Koustav Das Chakladar, Aliah University*", unsafe_allow_html=True)

# -------------------------
# Main app logic
# -------------------------
if not st.session_state.authenticated:
    login_page()
else:
    # Sidebar: simple logout
    st.sidebar.write(f"**Logged in as:** {st.session_state.username}")
    if st.sidebar.button("Logout"):
        st.session_state.authenticated = False
        st.session_state.username = None
        st.session_state.user_choice = {k: None for k in CHOICE_ORDER}
        st.session_state.last_problem = None
        st.session_state.output_folder = None
        st.session_state.page = "home"
        st.experimental_rerun()

    st.markdown(f"**Welcome, {st.session_state.username}** | {breadcrumb_text()}")
    page = st.session_state.get("page", "home")
    if page == "home":
        home_page()
    elif page == "time_selection":
        time_selection_page()
    elif page == "body_selection":
        body_selection_page()
    elif page == "potential_selection":
        potential_selection_page()
    elif page == "perturbation_selection":
        perturbation_selection_page()
    elif page == "perturbation_type_selection":
        perturbation_type_selection_page()
    elif page == "problem_selection":
        problem_selection_page()
    elif page == "input_page":
        input_page()
    elif page == "plot_page":
        plot_page()
    elif page == "theory_plot_page":
        theory_plot_page()

