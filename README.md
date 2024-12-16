# MillPlus
The Mill+ application user interface is structured into multiple tabs, each dedicated to specific features.

1. Tool material and modal parameters tab: This tab allows the user to input material-specific parameters and modal properties such as stiffness, damping, and natural frequency for both x and y directions. This data is required for simulating the vibrations during the milling process.

	2. Milling operation tab:  Users can set operational parameters, such as spindle speed, depth of cut, and feed per tooth. Simulation settings, including the number of simulation periods. The RUN button initiates the milling simulation, visually tracking progress with a gauge.

	3. Stability charts tab: this tab provides stability diagrams displaying peak-to-peak force (in x, y and cross xy components), power consumption (P_c), mass removal rate (MRR), and surface roughness (R_a) as functions of spindle speed and axial depth of cut. These multivariable plots help users assess the optimal process settings by providing visualization of stability loss.

	4. Dynamic forces tab: the dynamic forces tab calculates and displays force components (F_x and F_y) over time. Additionally, Fast Fourier Transform (FFT) analysis is provided to identify frequencies associated with chatter and assess vibration patterns. Besides, this tab can be combined with the Optimization tool tab, to assess the goodness of the ‘Optimum point’ in terms of frequency spectrum content.
