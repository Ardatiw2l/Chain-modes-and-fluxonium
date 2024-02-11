# _interactive_ploting_functions.py
import matplotlib.pyplot as plt

class InteractiveScatterPlot_dispestion:
    def __init__(self, x, y, data_values):
        # Initialize the matplotlib figure and axis objects.
        self.fig, self.ax = plt.subplots()

        # Create a scatter plot with 'x' markers for the given x and y data.
        self.sc = self.ax.scatter(x, y, marker='x')

        # Store the additional data values to be displayed on hover.
        self.data_values = data_values

        # Set up an annotation object (invisible at start) for displaying data values on hover.
        # The annotation is positioned with an offset from the hovered point.
        self.annot = self.ax.annotate("", xy=(0,0), xytext=(20,20), textcoords="offset points",
                                      bbox=dict(boxstyle="round", fc="w"),
                                      arrowprops=dict(arrowstyle="->"))
        self.annot.set_visible(False)

        # Connect the hover event handling function to the figure.
        self.fig.canvas.mpl_connect("motion_notify_event", self.hover)

        # Display the plot.
        self.update_plot()

    def update_annot(self, ind):
        # Update the position and text of the annotation based on the index of the hovered point.
        pos = self.sc.get_offsets()[ind["ind"][0]]
        self.annot.xy = pos
        text = f"g = {self.data_values[ind['ind'][0]]:.3f} GHz"
        self.annot.set_text(text)
        self.annot.get_bbox_patch().set_alpha(0.4)

    def hover(self, event):
        # Event handling function for mouse movement. Updates the annotation visibility based on hover.
        if not (event.inaxes == self.ax):
            # Hide the annotation if the mouse is not over the axis.
            if self.annot.get_visible():
                self.annot.set_visible(False)
                self.fig.canvas.draw_idle()
            return

        # Check if the mouse is over a scatter plot point.
        cont, ind = self.sc.contains(event)
        if cont:
            # If so, update and display the annotation.
            self.update_annot(ind)
            if not self.annot.get_visible():
                self.annot.set_visible(True)
                self.fig.canvas.draw_idle()
        else:
            # If not, hide the annotation.
            if self.annot.get_visible():
                self.annot.set_visible(False)
                self.fig.canvas.draw_idle()

    def update_plot(self):
        # Render the plot. This method can be used to redraw the plot if needed.
        plt.show()

    def set_title(self, title):
        # Set the title of the axis and redraw the figure.
        self.ax.set_title(title)
        self.fig.canvas.draw_idle()

