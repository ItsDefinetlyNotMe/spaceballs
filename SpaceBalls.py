import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import scrolledtext
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from tkinter import filedialog
import copy
import random

COLORS = ['b','g','r','c','m','y','k','w','blue','green','red','cyan','magenta','yellow','black','white','orange','purple','pink']

class Movable():
    def __init__(self, velocity, position, color='r'):
        self.velocity = velocity
        self.position = position
        if color in COLORS:
            self.color = color
        else:
            self.color = 'r'

    def __str__(self):
        return f"position: {self.position}, velocity {self.velocity}"

class Planet(Movable):
    def __init__(self, velocity, position, color='r'):
        super().__init__(velocity, position, color)


class Ball(Movable):
    def __init__(self, velocity, position):
        super().__init__(velocity, position)
        self.color = 'b'

    
    def collision(self, planet:Planet):
        self.velocity = 2 * planet.velocity - self.velocity
        self.position = 2 * planet.position - self.position

class SpaceBall():
    
    def __init__(self, epsilon, arguments):
        self.ball = None
        self.original_ball = None
        self.planets = None
        self.sequences = None
        self.active_sequence = None
        self.history = None
        self.epsilon = epsilon
        self.arguments = arguments

    def add_scene(self, path):
        """Prepares objects described in scene file for simulation"""
        with open(path, 'r') as f:
            lines = f.readlines()
            ball, planets, sequences, _, _ = self.parser(lines)
        self.ball = ball
        self.original_ball = copy.deepcopy(ball)
        self.planets = planets
        self.sequences = sequences
        self.history = {}
        for seq in sequences:
            self.history[str(seq)] = [(ball.velocity, ball.position, 0)]

    def random_scene(self, planet_amount, sequence_length):
        self.ball = Ball(random.uniform(-30, 30), random.uniform(-30, 30))
        self.original_ball = copy.deepcopy(self.ball)
        self.planets = {}
        self.history = {}
        for i in range(planet_amount):
            self.planets[f"P{i}"] = Planet(random.uniform(-40, 40), random.uniform(-40, 40))
        
        s = self.BuildValidSequences([], dict.fromkeys(list(self.planets.keys()), 3), self.ball, 0, same=False)
        s = [es for es in s if len(es) ==  sequence_length]

        if len(s) == 0:
            return False
        self.sequences = [s[-1]]
        #print(s[-1])

        for seq in self.sequences:
            self.history[str(seq)] = [(self.ball.velocity, self.ball.position, 0)]
        return True

    def parser(self, scene):
        """Parses .scene files
        
        Returns:
        --------
        Ball
            Ball object described in scene.
        [Planet]
            List of Planet objects described in the scene
        [[string]]
            List of all sequences described in the scene
        [int]
            List of all lines in which a scene is specified
        [int]
            List of all lines which are comments
            """
        planets = {}
        ball = None
        sequences = []
        seq_line = []
        comments = []

        for idx, line in enumerate(scene):
            if line.startswith('#'):
                comments.append(idx)
                continue
            elif line.startswith(' ') or line == '':
                continue

            arguments = line.split()   
            
            if len(arguments) >= 2: 
                if arguments[0] == 'Sequence':
                    sequences.append(arguments[1].split(','))
                    seq_line.append(idx)
                    continue
                velocity = float(arguments[1])
                position = float(arguments[2])
                if arguments[0] == 'Ball':
                    ball = Ball(velocity, position)    
                elif arguments[0] == 'Planet':
                    name = arguments[3]
                    if len(arguments) == 5:
                        color = arguments[4]
                        planets[name] = Planet(velocity, position, color)
                    else:
                        planets[name] = Planet(velocity, position)

        return ball, planets, sequences, seq_line, comments

    def time_of_collision(self, m1:Movable, m2:Movable):
        """returns time of collision between two objects"""
        return -(m1.position-m2.position)/(m1.velocity-m2.velocity)

    def is_collision_possible(self, m1:Movable, m2:Movable)->bool:
        """Determins if a collision is possible at the end of the active sequence"""
        return self.time_of_collision(m1, m2) >= self.history[self.active_sequence][-1][2] + self.epsilon

    def collision(self, planet):
        """adjusts the stats of the ball according to the collision"""
        time = self.time_of_collision(self.ball, planet)
        self.ball.collision(planet)
        self.history[self.active_sequence].append((self.ball.velocity, self.ball.position, time))

    def play_sequence(self, seq):
        """Creates a history by simulating a sequence"""
        self.active_sequence = str(seq)
        self.ball = copy.deepcopy(self.original_ball)
        self.history[self.active_sequence] = [(self.ball.velocity, self.ball.position, 0)]
        for c, collision in enumerate(seq):
            collision_planet = self.planets.get(collision, None)
            if collision_planet is None:
                return False, c
            if self.is_collision_possible(self.ball, collision_planet):
                self.collision(collision_planet)
            else:
                return False, c
        return True, None

    def play_sequences(self)->bool:
        """
        Simulates all sequences sequentialy
 
        Returns:
        --------
        bool
            True when every sequence is possible, false otherwise.
        Tupel[int, int]
            Index of Sequence and Planet that are not possible, None otherwise.

        Example:
        --------
        >>> space_ball.PlaySequence()
        True, None
        """
            
        for s, sequence in enumerate(self.sequences):
            is_valid, c = self.play_sequence(sequence)
            if not is_valid:
                return False, (s, c)
        return True, None
                
    def get_valid_planets(self, ball:Ball):
        valid_planets = []
        for p in self.planets.keys():
            if self.is_collision_possible(ball, self.planets[p]):
                valid_planets.append(p)
        return valid_planets

    def plot_plot(self):
        """
        Plots all sequences

        Returns:
        --------
        figure
            plt.figure of all sequences 
        """

        def get_axis(seq_amount, plot_amount):
            if seq_amount == 1 and plot_amount == 1:
                return axes

            if self.arguments['r']:
                if plots == 1:
                    return axes[idx]
                elif self.sequences == 1:
                    return axes[counter]
                else:
                    return axes[idx, counter]

            if seq_amount == 1:
                return axes[counter]
            elif plots == 1:
                return axes[idx]
            else:
                return axes[counter, idx]
            

        plots = [self.arguments['p'], self.arguments['d'], self.arguments['s'], self.arguments['v']].count(True)
        if self.arguments['r']:
            fig, axes = plt.subplots(len(self.sequences), plots, figsize=(5*plots, 5*len(self.sequences)))
        else:
            fig, axes = plt.subplots(plots, len(self.sequences), figsize=(5*len(self.sequences), 5*plots))

        for idx, seq in enumerate(self.sequences):
            ac_seq = str(seq)
            if len(self.history[ac_seq]) > 0:
                his = self.history[ac_seq]

                velocity = [item[0] for item in his]
                timestamps = [item[2] for item in his]
                positions = [item[1] for item in his]

                pwt = [p + t*v for p, t, v in zip(positions, timestamps, velocity)]
                last_time = timestamps[-1]
                planet_time = [0, last_time]

                counter = 0
                #what if len(sequence) and plots == 1? / what if len(seq)==0
                if self.arguments['p']:
                    ax1 = get_axis(len(self.sequences), plots)
                    
                    counter += 1
                    for x in self.arguments['hx']:
                        ax1.axvline(x=x, color='green', linestyle=':', linewidth=2)
                    #Plot 1
                    ax1.plot(timestamps, pwt, linestyle='-', color='b')
                    

                    for planet in self.planets.keys():
                        position = self.planets[planet].position
                        vel = self.planets[planet].velocity
                        pos = [position, position + vel*last_time]
                        ax1.plot(planet_time, pos, linestyle='-', color=self.planets[planet].color)

                    ax1.set_xlabel('Time (seconds)')
                    ax1.set_ylabel('Position')
                    ax1.set_title(f'Sequence {idx+1}: Position vs Time')
                    ax1.grid(True)
                    ax1.annotate(r'$t_{\text{end}}=$' + f'{round(timestamps[-1],3)}', xy=(0.95, 0.05), xycoords='axes fraction',ha='right',va='bottom')
                    ax1.annotate(r'$p_{\text{end}}=$' + f'{round(pwt[-1],3)}', xy=(0.95, 0.00), xycoords='axes fraction',ha='right',va='bottom')

                    if self.arguments['my'] is not None:
                        ax1.set_ylim(top=self.arguments['my'])
                    if self.arguments['ly'] is not None:
                        ax1.set_ylim(bottom=self.arguments['ly'])
                        
                if self.arguments['d']:
                    ax2 = get_axis(len(self.sequences), plots)
                    
                    for x in self.arguments['hx']:
                        ax2.axvline(x=x, color='green', linestyle=':', linewidth=2)
                    
                    counter += 1
                    #Plot 2
                    cummulative_position = []
                    prepos = pwt[0]
                    for i,p in enumerate(pwt):
                        distance = abs(p - prepos)
                        prepos = p
                        if i == 0:
                            cummulative_position.append(distance)
                            continue
                        cummulative_position.append(cummulative_position[-1] + distance)
                    ax2.plot(timestamps, cummulative_position, linestyle='-', color='b')

                    ax2.set_xlabel('Time (seconds)')
                    ax2.set_ylabel('Distance')
                    ax2.set_title(f'Sequence {idx+1}: Distance vs Time')
                    ax2.grid(True)
                    ax2.annotate(r'$s_{\text{end}}=$' + f'{round(cummulative_position[-1],3)}', xy=(0.95, 0.05), xycoords='axes fraction',ha='right',va='bottom')

                if self.arguments['v']:
                    ax3 = get_axis(len(self.sequences), plots)
  
                    for x in self.arguments['hx']:
                        ax3.axvline(x=x, color='green', linestyle=':', linewidth=2)
  
                    counter += 1
                    #Plot 3
                    ax3.step(timestamps, velocity, where='post', linestyle='-', color='b')

                    ax3.set_xlabel('Time (seconds)')
                    ax3.set_ylabel('Velocity')
                    ax3.set_title(f'Sequence {idx+1}: velocity vs Time')
                    ax3.grid(True)
                    ax3.annotate(r'$v_{\text{end}}=$' + f'{round(velocity[-1],3)}', xy=(0.95, 0.05), xycoords='axes fraction',ha='right',va='bottom')
                
                if self.arguments['s']:
                    ax4 = get_axis(len(self.sequences), plots)

                    for x in self.arguments['hx']:
                        ax4.axvline(x=x, color='green', linestyle=':', linewidth=2)

                    counter += 1
                    #Plot 3
                    ax4.step(timestamps, [abs(v) for v in velocity], where='post', linestyle='-', color='b')

                    ax4.set_xlabel('Time (seconds)')
                    ax4.set_ylabel('Speed')
                    ax4.set_title(f'Sequence {idx+1}: Speed vs Time')
                    ax4.grid(True)
                    ax4.annotate(r'$s_{\text{end}}=$' + f'{round(velocity[-1],3)}', xy=(0.95, 0.05), xycoords='axes fraction',ha='right',va='bottom')

        fig.tight_layout()
        return fig

    def show(self):
        """Shows current Plot"""
        plt.show()

    def print_details(self):
        """prints History of all sequences in terminal"""
        for key, history in self.history.items():
            print(f"Sequence: {key}")
            seq = key[1:-1:].split(',')
            h1 = history[0]
            history = history[1::]
            print(f"Initial condition: {h1[0]}*t + {h1[1]}")
            for col, ball_information in zip(seq,history):
                b_v, b_p, b_t = ball_information
                print(f"Collision with: {col}:\n\ttime: {b_t}\n\tball representation: {b_v}*t + {b_p}\n\tposition of collision: {b_v*b_t + b_p}")
            print("\n")
        print("------------------------------------")

    def BuildValidSequences(self, sequence, planets, ball, time, same=True):
        L = []
        L.append(sequence)
        for planetname in planets.keys():
            planet = self.planets[planetname]
            t=-1
            if ball.velocity != planet.velocity: 
                t = self.time_of_collision(ball, planet)
            if t >= time:
                s_prime = sequence + [planetname]
                if len(s_prime) > 1 and not same:
                    if s_prime[-1] == s_prime[-2]:
                        continue
                planets_prime = copy.deepcopy(planets)
                planets_prime[planetname] = planets_prime[planetname] - 1
                if planets_prime[planetname] <= 0:
                    planets_prime.pop(planetname, None)

                ball_prime =  Ball(2*planet.velocity - ball.velocity, 2*planet.position - ball.position)
            
                L += self.BuildValidSequences(s_prime, planets_prime, ball_prime, t, same)
        return L

    def newAlgorithm(self):
        seq = self.sequences[0]
        self.play_sequence(seq)
        time = self.history[str(seq)][-1][2]

        L = [seq]

        for n in range(len(seq)+1):
            for i in range(len(seq)+1):
                for j in range(i,len(seq)+1):
                    if n >= i and n <= j:
                        continue
                    else:
                        if n > j:
                            s = seq[:i] + seq[j:n] + seq[i:j] + seq[n::]
                            
                        elif n < i:
                            s = seq[:n] + seq[i:j] + seq[n:i] + seq[j::]


                        t,_ = self.play_sequence(s)
                        if t:
                            if round(self.history[str(s)][-1][2],3)==round(time,3):
                                L.append(s)
        return L

    def canShorten(self, sequences):
        for seq in sequences:
            for i in range(1, len(seq)):
                if seq[i] == seq[i-1]:
                    return True, seq
        return False, None

    def comparison(self):
        seq = self.sequences[0]
        _, bs = self.BuildSubShift()
        nA = self.newAlgorithm()
        c_nA , wor= self.canShorten(nA)
        print(c_nA)
        print(wor)
        nB = []
        for i in nA:
            if i not in nB:
                nB.append(i)
        nA = nB

        print(nA)
        exit()
        
        c_nA , wor= self.canShorten(nA)
        c_bs, wor2 = self.canShorten(bs)

        if not (c_nA and c_bs):
            if c_nA or c_bs:
                if not c_nA:
                    print("na")
                    output = f"#new Algorithm fails: \nBall {self.original_ball.velocity} {self.original_ball.position}\n"
                    output += "#Planets:\n"
                    for p in self.planets:
                        output += f"Planet {self.planets[p].velocity} {self.planets[p].position} {p}\n"
                    output += f"#Normal\nSequence {','.join(self.sequences[0])} \n#NF\nSequence {','.join(wor2)}\n-------------\n"
                    with open("NAF.txt", "a") as f:
                        f.write(output)
                elif not c_bs:
                    print("bs")
                    # Build output for file 2 (Subshift Algorithm fails)
                    output = f"#Subshift Algorithm fails: \nBall {self.original_ball.velocity} {self.original_ball.position}\n"
                    output += "#Planets:\n"
                    for p in self.planets:
                        output += f"Planet {self.planets[p].velocity} {self.planets[p].position} {p}\n"
                    output += f"#Normal\nSequence {','.join(self.sequences[0])} \n#NF\nSequence {','.join(wor)}\n-------------\n"
                    with open("SSF.txt", "a") as f:
                        f.write(output)
        
    def BuildSubShift(self):
        sequence = self.sequences[0]
        ball = self.original_ball

        #print(sequence)
        valid, _ = self.play_sequence(sequence)
        if not valid:
            return [], []
        else:
            vel = self.history[str(sequence)][-1][0]
            pos = self.history[str(sequence)][-1][1]
            time = self.history[str(sequence)][-1][2]

        R = set()
        samelen = []
        smallest = sequence
        SN = []
        planets = {} 
        for s in sequence:
            planets[s] = planets.get(s, 0) + 1
        
        t = 0
        L = self.BuildValidSequences(SN, planets, ball, t)
       
        #print(L)
        #print(len([l for l in L if len(l) == len(sequence)]))
        #print(len(L))
        L2 = []
        for seq in L:
            valid, _ = self.play_sequence(seq)
            if not valid:
                print("Error")
                exit()
            else:
                his = self.history[str(seq)][-1]
                nvel = his[0]
                npos = his[1]
                ntime = his[2]
            if round(nvel,3) == round(vel,3) and round(npos,3) == round(pos,3) and round(time,3)==round(ntime,3):
                L2.append(seq)
                if len(sequence) == len(seq):
                    samelen.append(seq)
                R.add(str(seq))
                if len(seq) < len(smallest):
                    smallest = seq

        #print("Subshift: " + str(R))
        print("same len: " + str(samelen))
        print("NF:" + str(smallest))

        return R, L2

    def form_normalform(self):
        sequence = self.sequences[0]
        planets = dict.fromkeys(set(sequence),0)
        self.ball = self.original_ball
        for idx, coll in enumerate(sequence):
            planets[coll] += -1 if idx % 2 == 0 else 1
        
        print(self._rec_normalform(seq=[],planets=planets, time=0))

    def _rec_normalform(self, seq, planets, time):
        ##calculate velocity,position of seq:
        total_velocity = 0
        total_position = 0

        for idx,col in enumerate(seq):
            idx = idx + 1
            if idx%2 == 0:
                total_velocity += 2*self.planets[col].velocity
                total_position += 2*self.planets[col].position
            else:
                total_velocity -= 2*self.planets[col].velocity
                total_position -= 2*self.planets[col].position
            print(f"TV: {total_velocity}")
        
        total_velocity += self.ball.velocity
        total_position += self.ball.position

        
        if len(seq) % 2 == 1:
            total_velocity *= -1
            total_position *= -1
        

        #find possible continuations
        pos_seqs = []
        b = 0
        for pl, amount in planets.items():
            planet = self.planets[pl]
            if abs(amount) > 0:
                b = 1
                #is it valid ?
                print(f"{total_velocity}*t + {total_position}")
                if (total_velocity-planet.velocity) == 0:
                    n_time = -1
                else:
                    n_time = (planet.position-total_position)/(total_velocity-planet.velocity)
                print(f"Seq: {seq}, ct: {time+self.epsilon}, {pl}: time: {n_time}")
                if  n_time >= time + self.epsilon:
                    n_planets = copy.deepcopy(planets)
                    n_planets[pl] = abs(amount) - 1
                    n_seq = copy.deepcopy(seq)
                    n_seq.append(pl)
                    pos_seqs.append(self._rec_normalform(n_seq, n_planets, n_time))
        
        if b==0:
            return seq

        return pos_seqs             

    def update_plot(self, canvas, fig):
        """Updates the plot in the UI"""
        for widget in canvas.winfo_children():
            widget.destroy()
        if fig:
            canvas_plot = FigureCanvasTkAgg(fig, master=canvas)
            canvas_plot.draw()
            canvas_plot.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

    def add_live_scene(self, path=None):
        """Tkinter interface with real-time plot update"""
        root = tk.Tk()
        root.title("SpaceBall Scene Editor")

        root.bind("<Destroy>", lambda event: plt.close(self.figure))

        text_frame = tk.Frame(root)
        text_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        text_box = scrolledtext.ScrolledText(text_frame, wrap=tk.WORD, width=50)
        text_box.pack(fill=tk.BOTH, expand=True)
        text_box.tag_configure("error", foreground="red")
        text_box.tag_configure("comment", foreground="green")

        right_frame = tk.Frame(root)
        right_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        valid_planet_frame = tk.Frame(right_frame)
        valid_planet_frame.pack(side=tk.BOTTOM, fill=tk.X)

        plot_frame = tk.Frame(right_frame)
        plot_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        valid_planet_label = tk.Label(valid_planet_frame, text="Valid Planets: ", anchor="w")
        valid_planet_label.pack(side=tk.LEFT, fill=tk.X, padx=10, pady=5)

        
        self.figure = None
     
        def on_text_change(event=None):
            cursor_position = text_box.index(tk.INSERT)
            scene_content = text_box.get("1.0", tk.END).strip()
            try:
                ball, planets, sequences, seq_lines, comments = self.parser(scene_content.split('\n'))

                self.ball = ball
                self.original_ball = copy.deepcopy(self.ball)
                self.planets = planets
                
                text_box.tag_remove("comment", "1.0", tk.END)
                for comment in comments:
                    text_box.tag_add("comment", f"{comment+1}.{0}", f"{comment+2}.{0}")
                    
                if sequences != [] and ball is not None:
                    self.sequences = sequences
                    self.history = {}
                    for seq in sequences:
                        self.history[str(seq)] = [(ball.velocity, ball.position, 0)]

                    text_box.tag_remove("error", "1.0", tk.END)
                    is_sim_complete, index = self.play_sequences()

                    if not is_sim_complete:
                        #calculate line index by sorting through all lines containing sequences and calculate their length
                        lines = scene_content.split('\n')
                        sequence_index = []
                        for i,l in enumerate(lines):
                            if l.startswith("Sequence"):
                                sequence_index.append(i)
                        err_line = sequence_index[index[0]]+1

                        err_char = 10 + len(",".join(self.sequences[index[0]][:index[1]]))
                        name_len = len(self.sequences[index[0]][index[1]])

                        text_box.tag_add("error", f"{err_line}.{0}", f"{err_line}.{8}")
                        text_box.tag_add("error", f"{err_line}.{err_char}", f"{err_line}.{err_char+name_len}")
                        valid_planet_label.config(text="Not a valid Sequence or Sequence was repeated")
                    else:
                        cursor_line = int(float(cursor_position)) - 1

                        if cursor_line in seq_lines:
                            s = self.sequences[seq_lines.index(cursor_line)]
                            self.active_sequence = str(s)
                            lt_v, lt_p, _ = self.history[str(s)][-1] 
                            valid_planets = self.get_valid_planets(Ball(lt_v, lt_p))
                            valid_planet_label.config(text="Sequence: " + str(seq_lines.index(cursor_line) + 1) + " Valid Planets: " + ", ".join(valid_planets))
                        else:
                            valid_planet_label.config(text="")
               
            except Exception as e:
                valid_planet_label.config(text=f"Error: {str(e)}")

        def ui_plot():
            plt.close(self.figure)    
            self.figure = self.plot_plot()
            self.update_plot(plot_frame, self.figure)


        if path is not None:
            with open(path,'r') as f:
                text_box.insert(tk.END, "".join(f.readlines()))
            on_text_change()
            ui_plot()


        text_box.bind("<KeyRelease>", on_text_change)

        def save_to_file():
            file_path = filedialog.asksaveasfilename(defaultextension=".scene", 
                                                    filetypes=[("Scene Files", "*.scene"), ("All Files", "*.*")])
            if file_path:
                content = text_box.get("1.0", tk.END).strip()
                with open(file_path, 'w') as file:
                    file.write(content)

        plot_button = tk.Button(valid_planet_frame, text="   Plot    ", command=ui_plot)
        plot_button.pack(side=tk.RIGHT, fill=tk.X, padx=10, pady=5)

        plot_button = tk.Button(valid_planet_frame, text="show plot in plt", command=self.show)
        plot_button.pack(side=tk.RIGHT, fill=tk.X, padx=10, pady=5)

        nf_button = tk.Button(valid_planet_frame, text="NF", command=self.comparison)#self.BuildSubShift)
        nf_button.pack(side=tk.RIGHT, fill=tk.X, padx=10, pady=5)

        print_button = tk.Button(valid_planet_frame, text="Print Details", command=self.print_details)
        print_button.pack(side=tk.RIGHT, fill=tk.X, padx=10, pady=5)

        save_button = tk.Button(text_frame, text="Save", command=save_to_file)
        save_button.pack(fill=tk.X, padx=5, pady=5)

        def on_toggle():
            self.arguments['s'] = speed_check.get()
            self.arguments['v'] = velocity_check.get()
            self.arguments['d'] = distance_check.get()
            self.arguments['r'] = reverse_check.get()

        speed_check = tk.BooleanVar(value=self.arguments["s"]) 
        velocity_check = tk.BooleanVar(value=self.arguments["v"]) 
        distance_check = tk.BooleanVar(value=self.arguments["d"]) 
        reverse_check = tk.BooleanVar(value=self.arguments["r"]) 

        s_check = tk.Checkbutton(text_frame, text="Plot Speed", variable=speed_check, command=on_toggle)
        s_check.pack(padx=20, pady=20, side=tk.LEFT)
        v_check = tk.Checkbutton(text_frame, text="Plot Velocity", variable=velocity_check, command=on_toggle)
        v_check.pack(padx=20, pady=20, side=tk.LEFT)
        d_check = tk.Checkbutton(text_frame, text="Plot Distance", variable=distance_check, command=on_toggle)
        d_check.pack(padx=20, pady=20, side=tk.LEFT)
        r_check = tk.Checkbutton(text_frame, text="Reverse Plot", variable=reverse_check, command=on_toggle)
        r_check.pack(padx=20, pady=20, side=tk.LEFT)

        root.mainloop()