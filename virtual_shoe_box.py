import tkinter as tk
from tkinter import ttk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# Settings

font_normal = ("Arial", 10)
font_bold = ("Arial", 10, "bold")
spinbox_width = 10

# Root window

root = tk.Tk()
#root.geometry('300x200')
root.resizable(False, False)
root.title('Virtual Shoe Box')

# Grid

for k in range(9):
    root.columnconfigure(k, weight=1)
for k in range(11):
    root.rowconfigure(k, weight=1)

# Source

lbl_src = tk.Label(root, text = "Source", font=font_bold)
lbl_src.grid(column=0, row=0, sticky = "nw", pady = 5, padx = 5)

lbl_src_d0 = tk.Label(root, text = "Distance [m]", font=font_normal)
lbl_src_d0.grid(column=1, row=0, sticky = "nw", pady = 5, padx = 5)

spinbox_src_d0_current_value = tk.StringVar(value=0)
spinbox_src_d0_current_value.set("1.0")
spinbox_src_d0 = ttk.Spinbox(
    root,
    from_=1,
    to=10,
    increment=0.1,
    textvariable=spinbox_src_d0_current_value,
    wrap=False,
    font=font_normal,
    width=spinbox_width)
spinbox_src_d0.grid(column=2, row=0, sticky = "new", pady = 5, padx = 5)

lbl_src_phi0 = tk.Label(root, text = "Azimuth [deg]", font=font_normal)
lbl_src_phi0.grid(column=3, row=0, sticky = "nw", pady = 5, padx = 5)

spinbox_src_phi0_current_value = tk.StringVar(value=0)
spinbox_src_phi0_current_value.set("0")
spinbox_src_phi0 = ttk.Spinbox(
    root,
    from_=0,
    to=359,
    increment=1,
    textvariable=spinbox_src_phi0_current_value,
    wrap=True,
    font=font_normal,
    width=spinbox_width)
spinbox_src_phi0.grid(column=4, row=0, sticky = "new", pady = 5, padx = 5)

lbl_src_theta0 = tk.Label(root, text = "Elevation [deg]", font=font_normal)
lbl_src_theta0.grid(column=5, row=0, sticky = "nw", pady = 5, padx = 5)

spinbox_src_theta0_current_value = tk.StringVar(value=0)
spinbox_src_theta0_current_value.set("0")
spinbox_src_theta0 = ttk.Spinbox(
    root,
    from_=-180,
    to=180,
    increment=1,
    textvariable=spinbox_src_theta0_current_value,
    wrap=True,
    font=font_normal,
    width=spinbox_width)
spinbox_src_theta0.grid(column=6, row=0, sticky = "new", pady = 5, padx = 5)

# Walls

lbl_walls = tk.Label(root, text = "Walls", font=font_bold)
lbl_walls.grid(column=0, row=1, sticky = "nw", pady = 5, padx = 5)

lbl_walll = tk.Label(root, text = "Left", font=font_normal)
lbl_walll.grid(column=0, row=2, sticky = "nw", pady = 5, padx = 5)

lbl_wallr = tk.Label(root, text = "Right", font=font_normal)
lbl_wallr.grid(column=0, row=3, sticky = "nw", pady = 5, padx = 5)

lbl_wallf = tk.Label(root, text = "Front", font=font_normal)
lbl_wallf.grid(column=0, row=4, sticky = "nw", pady = 5, padx = 5)

lbl_wallb = tk.Label(root, text = "Back", font=font_normal)
lbl_wallb.grid(column=0, row=5, sticky = "nw", pady = 5, padx = 5)

lbl_wallc = tk.Label(root, text = "Ceiling", font=font_normal)
lbl_wallc.grid(column=0, row=6, sticky = "nw", pady = 5, padx = 5)

lbl_wallo = tk.Label(root, text = "Floor", font=font_normal)
lbl_wallo.grid(column=0, row=7, sticky = "nw", pady = 5, padx = 5)

#

lbl_walls_d = tk.Label(root, text = "Distance [m]", font=font_normal)
lbl_walls_d.grid(column=1, row=1, sticky = "nw", pady = 5, padx = 5)

spinbox_walls_d_current_value = []
spinbox_walls_d = []
for k in range(6):
    spinbox_walls_d_current_value.append(tk.StringVar(value=0))
    spinbox_walls_d_current_value[k].set("2.0")
    spinbox_walls_d.append(ttk.Spinbox(
        root,
        from_=1,
        to=100,
        increment=0.1,
        textvariable=spinbox_walls_d_current_value[k],
        wrap=False,
        font=font_normal,
        width=spinbox_width))
    spinbox_walls_d[k].grid(column=1, row=k+2, sticky = "new", pady = 5, padx = 5)

#

lbl_walls_g = tk.Label(root, text = "Gain [linear]", font=font_normal)
lbl_walls_g.grid(column=2, row=1, sticky = "nw", pady = 5, padx = 5)

spinbox_walls_g_current_value = []
spinbox_walls_g = []
for k in range(6):
    spinbox_walls_g_current_value.append(tk.StringVar(value=0))
    spinbox_walls_g_current_value[k].set("1.0")
    spinbox_walls_g.append(ttk.Spinbox(
        root,
        from_=0,
        to=1,
        increment=0.1,
        textvariable=spinbox_walls_g_current_value[k],
        wrap=False,
        font=font_normal,
        width=spinbox_width))
    spinbox_walls_g[k].grid(column=2, row=k+2, sticky = "new", pady = 5, padx = 5)

#

lbl_walls_lf = tk.Label(root, text = "Low freq [Hz]", font=font_normal)
lbl_walls_lf.grid(column=3, row=1, sticky = "nw", pady = 5, padx = 5)

spinbox_walls_lf_current_value = []
spinbox_walls_lf = []
for k in range(6):
    spinbox_walls_lf_current_value.append(tk.StringVar(value=0))
    spinbox_walls_lf_current_value[k].set("100")
    spinbox_walls_lf.append(ttk.Spinbox(
        root,
        from_=10,
        to=1000,
        increment=1,
        textvariable=spinbox_walls_lf_current_value[k],
        wrap=False,
        font=font_normal,
        width=spinbox_width))
    spinbox_walls_lf[k].grid(column=3, row=k+2, sticky = "new", pady = 5, padx = 5)

#

lbl_walls_lg = tk.Label(root, text = "Low gain [dB]", font=font_normal)
lbl_walls_lg.grid(column=4, row=1, sticky = "nw", pady = 5, padx = 5)

spinbox_walls_lg_current_value = []
spinbox_walls_lg = []
for k in range(6):
    spinbox_walls_lg_current_value.append(tk.StringVar(value=0))
    spinbox_walls_lg_current_value[k].set("0.0")
    spinbox_walls_lg.append(ttk.Spinbox(
        root,
        from_=-20,
        to=20,
        increment=0.5,
        textvariable=spinbox_walls_lg_current_value[k],
        wrap=False,
        font=font_normal,
        width=spinbox_width))
    spinbox_walls_lg[k].grid(column=4, row=k+2, sticky = "new", pady = 5, padx = 5)

#

lbl_walls_lQ = tk.Label(root, text = "Low Q", font=font_normal)
lbl_walls_lQ.grid(column=5, row=1, sticky = "nw", pady = 5, padx = 5)

spinbox_walls_lQ_current_value = []
spinbox_walls_lQ = []
for k in range(6):
    spinbox_walls_lQ_current_value.append(tk.StringVar(value=0))
    spinbox_walls_lQ_current_value[k].set("0.71")
    spinbox_walls_lQ.append(ttk.Spinbox(
        root,
        from_=0,
        to=1,
        increment=0.01,
        textvariable=spinbox_walls_lQ_current_value[k],
        wrap=False,
        font=font_normal,
        width=spinbox_width))
    spinbox_walls_lQ[k].grid(column=5, row=k+2, sticky = "new", pady = 5, padx = 5)

#

lbl_walls_hf = tk.Label(root, text = "High freq [Hz]", font=font_normal)
lbl_walls_hf.grid(column=6, row=1, sticky = "nw", pady = 5, padx = 5)

spinbox_walls_hf_current_value = []
spinbox_walls_hf = []
for k in range(6):
    spinbox_walls_hf_current_value.append(tk.StringVar(value=0))
    spinbox_walls_hf_current_value[k].set("5000")
    spinbox_walls_hf.append(ttk.Spinbox(
        root,
        from_=100,
        to=20000,
        increment=10,
        textvariable=spinbox_walls_hf_current_value[k],
        wrap=False,
        font=font_normal,
        width=spinbox_width))
    spinbox_walls_hf[k].grid(column=6, row=k+2, sticky = "new", pady = 5, padx = 5)

#

lbl_walls_hg = tk.Label(root, text = "High gain [dB]", font=font_normal)
lbl_walls_hg.grid(column=7, row=1, sticky = "nw", pady = 5, padx = 5)

spinbox_walls_hg_current_value = []
spinbox_walls_hg = []
for k in range(6):
    spinbox_walls_hg_current_value.append(tk.StringVar(value=0))
    spinbox_walls_hg_current_value[k].set("0.0")
    spinbox_walls_hg.append(ttk.Spinbox(
        root,
        from_=-20,
        to=20,
        increment=0.5,
        textvariable=spinbox_walls_hg_current_value[k],
        wrap=False,
        font=font_normal,
        width=spinbox_width))
    spinbox_walls_hg[k].grid(column=7, row=k+2, sticky = "new", pady = 5, padx = 5)

#

lbl_walls_hQ = tk.Label(root, text = "High Q", font=font_normal)
lbl_walls_hQ.grid(column=8, row=1, sticky = "nw", pady = 5, padx = 5)

spinbox_walls_hQ_current_value = []
spinbox_walls_hQ = []
for k in range(6):
    spinbox_walls_hQ_current_value.append(tk.StringVar(value=0))
    spinbox_walls_hQ_current_value[k].set("0.71")
    spinbox_walls_hQ.append(ttk.Spinbox(
        root,
        from_=0,
        to=1,
        increment=0.01,
        textvariable=spinbox_walls_hQ_current_value[k],
        wrap=False,
        font=font_normal,
        width=spinbox_width))
    spinbox_walls_hQ[k].grid(column=8, row=k+2, sticky = "new", pady = 5, padx = 5)

# Air

lbl_air = tk.Label(root, text = "Air", font=font_bold)
lbl_air.grid(column=0, row=8, sticky = "nw", pady = 5, padx = 5)

lbl_air_abs = tk.Label(root, text = "Air absorption", font=font_normal)
lbl_air_abs.grid(column=1, row=8, sticky = "nw", pady = 5, padx = 5)

spinbox_air_abs_current_value = tk.StringVar(value=0)
spinbox_air_abs_current_value.set("On")
spinbox_air_abs = ttk.Spinbox(
    root,
    textvariable=spinbox_air_abs_current_value,
    wrap=True,
    font=font_normal,
    values=("On", "Off"),
    width=spinbox_width)
spinbox_air_abs.grid(column=2, row=8, sticky = "new", pady = 5, padx = 5)

lbl_air_h = tk.Label(root, text = "Humidity", font=font_normal)
lbl_air_h.grid(column=3, row=8, sticky = "nw", pady = 5, padx = 5)

spinbox_air_h_current_value = tk.StringVar(value=0)
spinbox_air_h_current_value.set("0.6")
spinbox_air_h = ttk.Spinbox(
    root,
    from_=0,
    to=1,
    increment=0.1,
    textvariable=spinbox_air_h_current_value,
    wrap=False,
    font=font_normal,
    width=spinbox_width)
spinbox_air_h.grid(column=4, row=8, sticky = "new", pady = 5, padx = 5)

lbl_air_v = tk.Label(root, text = "Speed [m/s]", font=font_normal)
lbl_air_v.grid(column=5, row=8, sticky = "nw", pady = 5, padx = 5)

spinbox_air_v_current_value = tk.StringVar(value=0)
spinbox_air_v_current_value.set("340")
spinbox_air_v = ttk.Spinbox(
    root,
    from_=300,
    to=400,
    increment=1,
    textvariable=spinbox_air_v_current_value,
    wrap=False,
    font=font_normal,
    width=spinbox_width)
spinbox_air_v.grid(column=6, row=8, sticky = "new", pady = 5, padx = 5)

# IR

lbl_IR = tk.Label(root, text = "IR", font=font_bold)
lbl_IR.grid(column=0, row=9, sticky = "nw", pady = 5, padx = 5)

lbl_IR_o = tk.Label(root, text = "Max order", font=font_normal)
lbl_IR_o.grid(column=1, row=9, sticky = "nw", pady = 5, padx = 5)

spinbox_IR_o_current_value = tk.StringVar(value=0)
spinbox_IR_o_current_value.set("2")
spinbox_IR_o = ttk.Spinbox(
    root,
    from_=0,
    to=100,
    increment=1,
    textvariable=spinbox_IR_o_current_value,
    wrap=False,
    font=font_normal,
    width=spinbox_width)
spinbox_IR_o.grid(column=2, row=9, sticky = "new", pady = 5, padx = 5)

lbl_IR_t = tk.Label(root, text = "Max time [s]", font=font_normal)
lbl_IR_t.grid(column=3, row=9, sticky = "nw", pady = 5, padx = 5)

spinbox_IR_t_current_value = tk.StringVar(value=0)
spinbox_IR_t_current_value.set("2.0")
spinbox_IR_t = ttk.Spinbox(
    root,
    from_=0,
    to=100,
    increment=0.1,
    textvariable=spinbox_IR_t_current_value,
    wrap=False,
    font=font_normal,
    width=spinbox_width)
spinbox_IR_t.grid(column=4, row=9, sticky = "new", pady = 5, padx = 5)

# Matplotlib

plotFrame = tk.Frame(root)
plotFrame.grid(column=0, row=10, columnspan=9, sticky = "ew", pady = 5, padx = 5)
f = Figure(figsize=(9, 2), dpi=100)

ax = f.add_subplot(121)
ax.set_xlim([-0.2, 1.2])
ax.set_ylim([-0.2, 1.2])

ax = f.add_subplot(122)
ax.set_xlim([-0.2, 1.2])
ax.set_ylim([-0.2, 1.2])

canvas = FigureCanvasTkAgg(f, master=plotFrame)
canvas.draw()
canvas.get_tk_widget().grid()

# Buttons

def IR_calc():
   pass

def IR_save():
   pass

def audio_load():
   pass

def audio_listen():
   pass

button_IR_calc = tk.Button(text ="Calculate IR", command = IR_calc)
button_IR_calc.grid(column=0, row=11, columnspan=2, sticky = "ew", pady = 5, padx = 5)

button_IR_save = tk.Button(text ="Save IR", command = IR_save)
button_IR_save.grid(column=2, row=11, columnspan=2, sticky = "ew", pady = 5, padx = 5)

button_audio_load = tk.Button(text ="Load audio file", command = audio_load)
button_audio_load.grid(column=4, row=11, columnspan=2, sticky = "ew", pady = 5, padx = 5)

button_audio_listen = tk.Button(text ="Listen to ...", command = audio_listen)
button_audio_listen.grid(column=6, row=11, columnspan=3, sticky = "ew", pady = 5, padx = 5)

# Main loop
root.mainloop()
