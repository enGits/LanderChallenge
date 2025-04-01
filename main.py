import arcade
import math
import sys
import datetime
import numpy as np
import random
import platform

SCREEN_WIDTH  = 0
SCREEN_HEIGHT = 0

IMAGE_PREFIX = ''

if hasattr(sys, '_MEIPASS'):
    IMAGE_PREFIX = sys._MEIPASS + '/'

def get_screen_size():
    system = platform.system()
    if system == "Windows":
        # Do NOT set DPI awareness here!
        try:
            import ctypes
            user32 = ctypes.windll.user32
            width = user32.GetSystemMetrics(0)  # SM_CXSCREEN
            height = user32.GetSystemMetrics(1) # SM_CYSCREEN
            return width, height
        except Exception as e:
            print(f"Failed to get scaled screen size via ctypes: {e}")
    else:
        # On other platforms, fallback to arcade or tkinter
        try:
            import tkinter as tk
            root = tk.Tk()
            root.withdraw()
            width = root.winfo_screenwidth()
            height = root.winfo_screenheight()
            root.destroy()
            return width, height
        except Exception as e:
            print(f"Fallback to arcade.get_display_size() due to: {e}")

    import arcade
    return arcade.get_display_size()

# Get screen dimensions
SCREEN_WIDTH, SCREEN_HEIGHT = get_screen_size()

# Get screen dimensions
SCREEN_WIDTH, SCREEN_HEIGHT = get_screen_size()

# Optional: Apply margin so the window always fits
MARGIN         = 70
SCREEN_WIDTH  -= MARGIN
SCREEN_HEIGHT -= MARGIN

SCREEN_TITLE   = 'enGits Lunar Lander Challenge'
NUM_MOUNTAINS  = 1000
ALT0           = 2.0
GAIN           = 5.0
ERROR_TIME     = 5.0
ZOOM_LENGTH    = 0.25
FEAT_THRESHOLD = 1000

HUD_FONT       = 'Courier New'
HUD_FONT_SIZE  = 10

MOON_COLOR     = (100,100,100)
MOUNTAIN_COLOR = (80,80,80)
GREEN          = (0, 200, 0)
RED            = arcade.color.LIGHT_RED_OCHRE
ORANGE         = arcade.color.ORANGE



class SurfaceFeature():
    
    def __init__(self, planet, theta=0):
        self.theta  = theta
        self.color  = arcade.color.RED
        self.planet = planet
        self.planet.feature_list.append(self)
        self.points = []
        
    def draw(self):
        alpha = self.planet.alpha + math.radians(self.theta)
        ref   = self.planet.game.reference
        x0    = ref.x
        y0    = ref.y
        pts   = []
        #
        x_min = -FEAT_THRESHOLD
        x_max =  FEAT_THRESHOLD + SCREEN_WIDTH
        y_min = -FEAT_THRESHOLD
        y_max =  FEAT_THRESHOLD + SCREEN_HEIGHT
        #
        draw = False
        for X in self.points:
            x1  = X[0]
            y1  = X[1] + self.planet.radius
            x2  = x1 * math.cos(alpha) - y1 * math.sin(alpha)
            y2  = x1 * math.sin(alpha) + y1 * math.cos(alpha)
            x2 += self.planet.x - x0
            y2 += self.planet.y - y0
            x3  = x2 * math.cos(ref.alpha) - y2 * math.sin(ref.alpha)
            y3  = x2 * math.sin(ref.alpha) + y2 * math.cos(ref.alpha)
            _x  = x3 * self.planet.game.scaleFactor() + SCREEN_WIDTH/2
            _y  = y3 * self.planet.game.scaleFactor() + SCREEN_HEIGHT/2
            if _x > x_min and _x < x_max and _y > y_min and _y < y_max:
                draw = True
            pts.append((_x, _y))
        if draw:
            arcade.draw_polygon_filled(pts, self.color)

    
class Mountain(SurfaceFeature):
    
    def __init__(self, planet, height, width, middle_width=None, theta=0):  
        super().__init__(planet, theta)
        self.color  = MOUNTAIN_COLOR
        if middle_width is None:
            middle_width = width/2
        self.points.append((       -width/2, -100))
        self.points.append((-middle_width/2,  height/2))
        self.points.append((              0,  height))
        self.points.append(( middle_width/2,  height/2))
        self.points.append((        width/2, -100))
        
        
        

class Body():
    
    G = 6.67430e-11  # Gravitational constant in m^3 kg^-1 s^-2
    
    def __init__(self, game):
        self.x     = 0.0  # Position in x-direction
        self.y     = 0.0  # Position in y-direction
        self.u     = 0.0  # Velocity component in x-direction
        self.v     = 0.0  # Velocity component in y-direction
        self.alpha = 0.0  # Angle of the body in radians
        self.omega = 0.0  # Rotational velocity (radiants per second)
        self.mass  = 0.0  # Mass of the body in kg
        self.game  = game
        self._x = 0.0
        self._y = 0.0
        self._a = 0.0
        self.acc_x        = 0.0
        self.acc_y        = 0.0
        self.acc_omega    = 0.0
        self.alpha        = 0.0
        self.color        = arcade.color.LIGHT_GRAY
        self.radius       = 0.0
        self.sprite       = None
        self.trajectory   = None
        self.scale_factor = 1.0
        self.name         = 'unnamed'
        self.docked_to    = None
        self.feature_list = []
        game.body_list.append(self)
        
    @staticmethod
    def safe_angle(dx, dy):
        """
        Compute angle in [0, 2π) from dx and dy.
        """
        angle = math.atan2(dy, dx)
        return angle
        
    def compute_orbit_2d(self, body, N = 101):
        G = Body.G
        m = body.mass
        r = [self.x - body.x, self.y - body.y]
        v = [self.u - body.u, self.v - body.v]        
        r = np.array(r, dtype=float)
        v = np.array(v, dtype=float)

        r_mag = np.linalg.norm(r)
        v_mag = np.linalg.norm(v)

        # Specific mechanical energy
        epsilon = 0.5 * v_mag**2 - G * m / r_mag

        # Angular momentum (scalar in 2D)
        h = r[0]*v[1] - r[1]*v[0]

        # Semi-major axis
        a = -G * m / (2 * epsilon)

        # Eccentricity vector
        e_vec = (1 / (G * m)) * np.array([
            (v_mag**2 - G * m / r_mag) * r[0] - np.dot(r, v) * v[0],
            (v_mag**2 - G * m / r_mag) * r[1] - np.dot(r, v) * v[1]
        ])
        e = np.linalg.norm(e_vec)

        # True anomaly (angle between position and periapsis)
        theta = np.arccos(np.dot(e_vec, r) / (e * r_mag))
        if np.dot(r, v) < 0:  # if moving away from periapsis
            theta = 2*np.pi - theta

        # Orbit equation as a function
        # def r_of_theta(th):
        #     return h**2 / (G * m * (1 + e * np.cos(th)))
        
        # Sample theta values over the orbit
        if e < 1.0:
            # Elliptical orbit: sample full 0..2pi
            theta_samples = np.linspace(0, 2*np.pi, N)
        else:
            # Hyperbolic orbit: sample a reasonable range around periapsis
            max_theta = np.arccos(-1 / e) if e > 1 else np.pi
            theta_samples = np.linspace(-max_theta, max_theta, N)

        # Orbit shape in polar coordinates
        r_samples = h**2 / (G * m * (1 + e * np.cos(theta_samples)))

        # Convert to Cartesian coordinates in orbital plane
        x_orb = r_samples * np.cos(theta_samples)
        y_orb = r_samples * np.sin(theta_samples)
        orbit_points = np.stack((x_orb, y_orb), axis=-1)

        # Rotate orbit so periapsis aligns with eccentricity vector
        angle_to_x = np.arctan2(e_vec[1], e_vec[0])
        rotation_matrix = np.array([
            [np.cos(angle_to_x), -np.sin(angle_to_x)],
            [np.sin(angle_to_x),  np.cos(angle_to_x)]
        ])
        self.trajectory = orbit_points @ rotation_matrix.T
        
    def draw_high_res_circle(self, x, y, radius, color, num_segments=120):
        #
        # get the minimum and maximum angles
        #
        min_angle = 0.0
        max_angle = 2*math.pi
        #
        # 4 ----- 3
        # |       |
        # 1 ----- 2
        #
        if y < 0 or y > SCREEN_HEIGHT or x < 0 or x > SCREEN_WIDTH:
            angle1 = Body.safe_angle(-x, -y)
            angle2 = Body.safe_angle(SCREEN_WIDTH - x, -y)
            angle3 = Body.safe_angle(SCREEN_WIDTH - x, SCREEN_HEIGHT - y)
            angle4 = Body.safe_angle(-x, SCREEN_HEIGHT - y)
            if x > SCREEN_WIDTH:
                if angle1 < 0:
                    angle1 += 2*math.pi
                if angle2 < 0:
                    angle2 += 2*math.pi
                if angle3 < 0:
                    angle3 += 2*math.pi
                if angle4 < 0:
                    angle4 += 2*math.pi
            min_angle = min(angle1, angle2, angle3, angle4)
            max_angle = max(angle1, angle2, angle3, angle4)
        #
        # draw the circle in segments
        #
        points = [(x, y)]  # Center of the circle
        for i in range(num_segments + 1):
            angle = min_angle + (max_angle - min_angle) * i / num_segments
            #angle = 2 * math.pi * i / num_segments
            px = x + radius * math.cos(angle)
            py = y + radius * math.sin(angle)
            points.append((px, py))
        arcade.draw_polygon_filled(points, color)        
        
    def update_metrics(self):
        if self == self.game.reference:
            self._x = SCREEN_WIDTH/2
            self._y = SCREEN_HEIGHT/2
            self._a = 0.0
            return
        x0 = self.game.reference.x
        y0 = self.game.reference.y
        a0 = self.game.reference.alpha
        #
        x1 = self.x - x0
        y1 = self.y - y0
        r1 = math.sqrt(x1**2 + y1**2)
        b1 = 0.0
        if r1 > 1:
            b1 = math.asin(abs(x1/r1))            
            if y1 < 0:
                b1 = math.pi - b1
            if x1 < 0:
                b1 = -b1
        b2 = b1 - a0
        x2 = r1 * math.sin(b2)
        y2 = r1 * math.cos(b2)
        self._a = math.degrees(self.alpha - a0)
        self._x = x2 * self.game.scaleFactor() + SCREEN_WIDTH/2
        self._y = y2 * self.game.scaleFactor() + SCREEN_HEIGHT/2
        
        #self._a = math.degrees(self.alpha)
        
    def update_physics(self, dt=0, Fx=0.0, Fy=0.0):
        for body in self.game.body_list:
            if body == self:
                continue
            if body.docked_to == self:
                continue
            if self.docked_to == body:
                continue
            dx = body.x - self.x
            dy = body.y - self.y
            r2 = dx**2 + dy**2
            r = math.sqrt(r2)
            if r2 == 0:
                continue
            r   = math.sqrt(r2)
            F   = (Body.G * self.mass * body.mass) / r2
            Fx += F * dx / r
            Fy += F * dy / r
        self.acc_x = Fx / self.mass
        self.acc_y = Fy / self.mass
        u = self.u
        v = self.v
        self.u += self.acc_x * dt
        self.v += self.acc_y * dt
        u = 0.5*(u + self.u)
        v = 0.5*(v + self.v)
        self.x += u * dt
        self.y += v * dt
        omega = self.omega
        self.omega += self.acc_omega * dt
        omega = 0.5*(omega + self.omega)
        self.alpha += omega * dt
        
    def draw(self):
        if self.game.dbg:
            pass
        if self.sprite is None:
            self.draw_high_res_circle(self._x, self._y, self.game.scaleFactor()*self.radius, self.color, num_segments=150)
            if self.game.reference != self.game.planet:
                for feature in self.feature_list:
                    feature.draw()
        else:
            pass

        
class Spacecraft(Body):
    def __init__(self, filename, scale, game):
        super().__init__(game)
        self.sprite = arcade.Sprite(filename, scale)
        for i in range(1,5):
            _file = filename.replace('01', '{:0>2d}'.format(i+1))
            self.sprite.append_texture(arcade.load_texture(IMAGE_PREFIX + _file))
        self.sprite.set_texture(0)
        game.sprite_list.append(self.sprite)
        self.thrust_level = 0
        self.ISP          = 314.0
        self.max_thrust   = 0.0
        self.fuel         = 0.0
        self.max_fuel     = 1.0
        self.dry_mass     = 0.0
        self.dock_dist    = 0.0
        self.target       = None
        self.landed       = False
        self.zoom_factor  = 1.0
        #
        self.tgt_vvel     = 0.0
        self.auto_pilot   = False
        self.errI         = 0.0
        self.crashed      = False

    def update_metrics(self):
        if self.game.reference != self.game.planet:
            d = 0.0
            if self.game.orbit_zoom:
                x1 = self.game.lander.x
                y1 = self.game.lander.y
                x2 = self.game.spacecraft.x
                y2 = self.game.spacecraft.y
                d  = max(1.0, math.sqrt((x1 - x2)**2 + (y1 - y2)**2))
                if d < 20.0 or self.game.lander.docked_to is not None:
                    self.scale_factor = self.game.real_zoom_factor
                else:
                    old_scale = self.scale_factor
                    new_scale = min(ZOOM_LENGTH*min(SCREEN_WIDTH, SCREEN_HEIGHT) / d, self.game.real_zoom_factor)
                    if new_scale < 0.5*old_scale:
                        new_scale = 0.5*old_scale
                    elif new_scale > 2.0*old_scale:
                        new_scale = 2.0*old_scale
                    else:
                        new_scale = old_scale
                    self.scale_factor = new_scale
            else:
                r   = math.sqrt(self.x**2 + self.y**2)
                R   = self.game.planet.radius
                alt = r - R
                self.scale_factor =  min(self.zoom_factor*ZOOM_LENGTH*min(SCREEN_WIDTH, SCREEN_HEIGHT) / alt, self.game.real_zoom_factor)
        super().update_metrics()
        if self.docked_to is not None:
            a  = math.radians(self.docked_to._a)
            dx = math.sin(a) * self.offset
            dy = math.cos(a) * self.offset
            self._x = self.docked_to._x + dx
            self._y = self.docked_to._y + dy
            self.alpha = self.docked_to.alpha + math.pi
        self.sprite.center_x = self._x
        self.sprite.center_y = self._y
        self.sprite.angle = self._a

    def draw(self):
        if self.game.dbg:
            pass
        #
        # draw HUD if required
        #
        if self == self.game.control_craft: # and self == self.game.reference:
            fs  = HUD_FONT_SIZE
            vel = math.sqrt(self.u**2 + self.v**2)
            alt = math.sqrt(self.x**2 + self.y**2) - self.game.planet.radius - ALT0
            x   = self._x  + 20
            y   = self._y
            fn  = HUD_FONT
            col = GREEN
            #
            if alt > 5.0e3:
                if self == self.game.lander:
                    col = ORANGE
                arcade.draw_text('km/h {:.2f}'.format(vel*3.6), x, y + fs + 2, col, fs, font_name=fn)
                arcade.draw_text('alt  {:.2f}'.format(alt/1000), x, y, col, fs, font_name=fn)
                arcade.draw_text('th.  {:.2f}%'.format(self.thrust_level), x, y - fs - 2, col, fs, font_name=fn)
                #arcade.draw_text('vv.  {:.2f}'.format(self.verticalVelocity()*3.6), x, y - 2*fs - 4, col, fs, font_name=fn)  
            else:
                if self == self.game.lander:
                    col = RED
                if vel < 5.0:
                    col = GREEN
                arcade.draw_text('m/s {:.2f}'.format(vel), x, y + fs + 2, col, fs, font_name=fn)
                arcade.draw_text('m   {:.2f}'.format(alt), x, y, col, fs, font_name=fn)
                arcade.draw_text('th. {:.2f}%'.format(self.thrust_level), x, y - fs - 2, col, fs, font_name=fn)
                #arcade.draw_text('vv.  {:.2f}'.format(self.verticalVelocity()), x, y - 2*fs - 4, col, fs, font_name=fn)
            #
            if self.auto_pilot:
                arcade.draw_text('tgt {:.2f}'.format(self.tgt_vvel), x, y - 3*fs - 6, col, fs, font_name=fn)
            u1 = self.u
            v1 = self.v
            u2 = u1
            v2 = v1
            u1_crit = 1.0
            u2_crit = u1_crit
            if self.target is not None:
                u2 -= self.target.u
                v2 -= self.target.v                
                u2_crit = 0.01
            u1_abs = math.sqrt(u1**2 + v1**2)
            u2_abs = math.sqrt(u2**2 + v2**2)
            if u1_abs > u1_crit:
                dx = u1 / u1_abs
                dy = v1 / u1_abs
                a  = Body.safe_angle(dx, dy)
                a += self.game.reference.alpha
                dx = 100*math.cos(a)
                dy = 100*math.sin(a)
                arcade.draw_line(self._x, self._y, self._x + dx, self._y + dy, col, 2)
            if u2_abs > u2_crit:
                dx = u2 / u2_abs
                dy = v2 / u2_abs
                a  = Body.safe_angle(dx, dy)
                a += self.game.reference.alpha
                dx = 100*math.cos(a)
                dy = 100*math.sin(a)
                arcade.draw_line(self._x, self._y, self._x + dx, self._y + dy, col, 1)
            
        #
        # draw orbital trajectory
        #
        if self.trajectory is not None and self.game.reference == self.game.planet:
            pts = []
            for x, y in self.trajectory:
                _x = (x - self.game.reference.x) * self.game.scaleFactor() + SCREEN_WIDTH/2
                _y = (y - self.game.reference.y) * self.game.scaleFactor() + SCREEN_HEIGHT/2
                pts.append((_x, _y))
            arcade.draw_line_strip(pts, GREEN, 1)
            n = int(len(pts)/2)
            arcade.draw_line(pts[0][0], pts[0][1], pts[n][0], pts[n][1], GREEN, 1)            
            km1 = (pts[0][0] - self.game.planet._x)**2 + (pts[0][1] - self.game.planet._y)**2
            km1 = 1e-3 * (math.sqrt(km1) / self.game.scaleFactor() - self.game.planet.radius)
            km2 = (pts[n][0] - self.game.planet._x)**2 + (pts[n][1] - self.game.planet._y)**2
            km2 = 1e-3 * (math.sqrt(km2) / self.game.scaleFactor() - self.game.planet.radius)
            arcade.draw_text('{:.1f}km'.format(km1), pts[0][0], pts[0][1], GREEN, 10)
            arcade.draw_text('{:.1f}km'.format(km2), pts[n][0], pts[n][1], GREEN, 10)
    
    def update_physics(self, dt=0):
        Fx = 0.0
        Fy = 0.0
        if self.name == 'Main Craft':
            pass
        update_trajectory = False
        if self.fuel > 0:
            F = self.max_thrust * 0.25*self.thrust_level
            self.fuel  = max(0.0, self.fuel - F / (9.81*self.ISP) * dt)
            Fx = F * math.sin(self.alpha)
            Fy = F * math.cos(self.alpha)
            if self.fuel < 1e-4:
                self.sprite.set_texture(0)
                self.thrust_level = 0
                self.fuel = 0.0
            if self.thrust_level > 0:
                update_trajectory = True
        alt = math.sqrt((self.x - self.game.planet.x)**2 + (self.y - self.game.planet.y)**2) - self.game.planet.radius
        if alt < ALT0:
            if not self.landed or self.thrust_level < 0.1:          
                u_abs = math.sqrt(self.u**2 + self.v**2)
                if u_abs > 5.0:
                    # CRASH if landing too fast
                    print("Crashed! Landing too fast! " + str(u_abs))
                    self.crashed = True
                self.u = 0
                self.v = 0
                self.acc_x = 0
                self.acc_y = 0
                self.omega = 0
                self.acc_omega = 0
                dx = self.x - self.game.planet.x
                dy = self.y - self.game.planet.y
                r   = math.sqrt(dx**2 + dy**2)
                dx /= r
                dy /= r
                beta  = math.degrees(math.atan2(dx, dy))
                alpha = math.degrees(self.alpha)
                if abs((alpha - beta + 180) % 360 - 180) > 20:
                    # CRASH due to bad landing angle
                    print("Crashed! Bad landing angle! " + str(abs((alpha - beta + 180) % 360 - 180)))
                    self.crashed = True
                self.x = self.game.planet.x + dx*(self.game.planet.radius + ALT0)
                self.y = self.game.planet.y + dy*(self.game.planet.radius + ALT0)
                self.trajectory = None
                self.thrust_level = 0
                self.sprite.set_texture(0)
                self.landed = True
                    
        if update_trajectory and self.trajectory is not None:
            self.compute_orbit_2d(self.game.planet)
        super().update_physics(dt, Fx, Fy)
        
    def resetAutoPilot(self):
        self.errI = 0.0
        self.tgt_hvel = 0.0        
        
    def verticalVelocity(self):
        dx  = self.x - self.game.planet.x
        dy  = self.y - self.game.planet.y
        r   = math.sqrt(dx**2 + dy**2)
        dx /= r
        dy /= r
        return self.u * dx + self.v * dy
        
    def altitude(self):
        dx  = self.x - self.game.planet.x
        dy  = self.y - self.game.planet.y
        r   = math.sqrt(dx**2 + dy**2)
        return r - self.game.planet.radius
        
    def autoPilot(self, dt):
        if self.auto_pilot:     
            if self.altitude() > 5000.0:
                self.auto_pilot = False
                return       
            err = self.tgt_vvel - self.verticalVelocity()
            self.errI += err * dt
            self.thrust_level = min(100.0, max(0.0, GAIN*(err + self.errI/ERROR_TIME)))                    
    
    def increase_thrust(self, amount):
        if not self.auto_pilot:            
            self.thrust_level = min(100, self.thrust_level + amount)
            idx = int(round(4*self.thrust_level/100))
            self.sprite.set_texture(idx)
        
    def decrease_thrust(self, amount):
        if not self.auto_pilot:
            self.thrust_level = max(0, self.thrust_level - amount)
            idx = int(round(4*self.thrust_level/100))
            self.sprite.set_texture(idx)
        
    def undock(self):
        if self.docked_to is not None:
            dx = math.sin(self.docked_to.alpha)
            dy = math.cos(self.docked_to.alpha)
            self.x = self.docked_to.x + dx*self.dock_dist
            self.y = self.docked_to.y + dy*self.dock_dist
            self.u = self.docked_to.u + 0.1*dx
            self.v = self.docked_to.v + 0.1*dy
            self.omega = self.docked_to.omega
            self.alpha = self.docked_to.alpha + math.pi
            self.acc_omega = self.docked_to.acc_omega
            self.acc_x = self.docked_to.acc_x
            self.acc_y = self.docked_to.acc_y
        self.docked_to = None
        
    def dock(self):
        csm = self.game.spacecraft
        dx  = self.x - csm.x
        dy  = self.y - csm.y
        a1  = math.degrees(self.alpha)
        a2  = math.degrees(csm.alpha)
        dist = math.sqrt(dx**2 + dy**2)
        if abs(dist - self.dock_dist) < 2.0:
            while a1 > 360:
                a1 -= 360
            while a1 < 0:
                a1 += 360
            while a2 > 360:
                a2 -= 360
            while a2 < 0:
                a2 += 360
            if abs(a1 - a2) - 180 < 10.0:
                self.docked_to = csm
                self.u = csm.u
                self.v = csm.v
                self.omega = csm.omega
                self.acc_omega = csm.acc_omega
                self.acc_x = csm.acc_x
                self.acc_y = csm.acc_y
                self.alpha = csm.alpha + math.pi
                self.x = csm.x + dx/dist*self.dock_dist
                self.y = csm.y + dy/dist*self.dock_dist
                self.game.control_craft = csm
                self.trajectory = None
    

class OrbitGame(arcade.View):
    def __init__(self, level=1):
        # super().__init__(SCREEN_WIDTH, SCREEN_HEIGHT, SCREEN_TITLE)
        super().__init__()

        self.time_factor      = 1.0
        self.sim_dt           = 0.1
        self.sim_time         = 0.0
        self.orbit_zoom       = False
        self.real_zoom_factor = 5.74

        arcade.set_background_color(arcade.color.BLACK)
        self.stars = [(random.randint(0, SCREEN_HEIGHT), random.randint(0, SCREEN_WIDTH)) for _ in range(100)]  # 100 stars
        self.body_list = []
        self.sprite_list = arcade.SpriteList()

        self.game_running = True 
        self.paused = False  
        self.game_over_view = None 
        self.show_crashed = False 

        # Planet parameters
        self.planet = Body(self)        
        self.planet.radius       = 1737.4e3      # Moon radius in metres
        self.planet.mass         = 7.34767309e22 # Moon mass in kg
        self.reference           = self.planet
        self.scale_factor        = 0.25*min(SCREEN_WIDTH, SCREEN_HEIGHT) / self.planet.radius
        # self.initial_scale       = self.scale
        self.planet.scale_factor = self.scale_factor
        self.planet.color        = MOON_COLOR
        self.planet.name         = 'Moon'
        self.dbg                 = False
        
        # Create random mountains on the Lunar surface
        # height H varies between 100m and 6000m
        # width W varies between 1*H and 3*H
        # middle width varies between 0.25*W and 0.75*W
        H_min = 100
        H_max = 6000
        for i in range(NUM_MOUNTAINS):
            H = H_min + (H_max - H_min)*np.random.rand()
            W = H + 2*H*np.random.rand()
            M = 0.25*W + 0.5*W*np.random.rand()
            theta = 360*np.random.rand()
            Mountain(self.planet, H, W, M, theta)
        

        # Spacecraft setup (roughly Apollo CSM for lunar mission)
        self.spacecraft = Spacecraft("craft_01.png", scale=0.1, game=self)
        self.spacecraft.dry_mass     = 8.3e3
        self.spacecraft.max_thrust   = 90e3
        self.spacecraft.fuel         = 11e3
        self.spacecraft.max_fuel     = self.spacecraft.fuel
        self.spacecraft.ISP          = 330.0
        self.spacecraft.y            = self.planet.radius + 300e3
        self.spacecraft.u            = math.sqrt(Body.G * self.planet.mass / self.spacecraft.y)
        self.spacecraft.scale_factor = 10*self.scale_factor
        self.spacecraft.name         = 'Main Craft'
        self.spacecraft.update_metrics()
        self.control_craft = self.spacecraft
               
        # Lander setup
        self.lander = Spacecraft("lander_01.png", scale=0.05, game=self)
        self.lander.dry_mass     = 3e3
        self.lander.max_thrust   = 30e3
        self.lander.fuel         = 11e3
        self.lander.max_fuel     = self.lander.fuel
        self.lander.ISP          = 330.0
        self.lander.y            = self.spacecraft.y
        self.lander.u            = self.spacecraft.u
        self.lander.alpha        = math.pi
        self.lander.docked_to    = self.spacecraft
        self.lander.dock_dist    = 5.6
        self.lander.offset       = 32.0
        self.lander.scale_factor = 20*self.scale_factor
        self.lander.name         = 'Lander'
        self.lander.update_metrics()
        
        arcade.schedule(self.on_update, 1 / 30)
        
        self.level = level
        if self.level == 2:
            self.set_level2()
        elif self.level == 3:
            self.set_level3()

    def set_level2(self):
        L = self.lander
        self.control_craft = L
        self.reference = L
        L.u = 0
        L.v = 0
        L.alpha = 0
        L.omega = 0
        L.acc_omega = 0
        L.acc_x     = 0
        L.acc_y     = 0
        L.x         = 0
        L.y         = -(self.planet.radius + 10e3)
        L.docked_to = None
        L.fuel     *= 0.5
        L.dock()

    def set_level3(self):
        L = self.lander
        self.control_craft = L
        self.reference = L
        L.y         = self.planet.radius + 100e3
        L.fuel     *= 0.1
        L.alpha     = 0
        L.omega     = 0
        L.acc_omega = 0
        L.acc_x     = 0
        L.acc_y     = 0
        L.x         = 0
        L.docked_to = None
        L.u         = math.sqrt(Body.G * self.planet.mass / L.y)

    def scaleFactor(self):
        return self.scale_factor
        
    def on_draw(self):
        self.clear()
        
        self.scale_factor = self.reference.scale_factor

        for star in self.stars:
            arcade.draw_circle_filled(star[0], star[1], 1, arcade.color.WHITE)  # Small white dots

        # update all metrics
        for body in self.body_list:
            body.update_metrics()
            
        # draw all non sprite bodies
        for body in self.body_list:
            if body.trajectory is None:
                body.draw()
        for body in self.body_list:
            if body.trajectory is not None:
                body.draw()

        # Draw spacecraft (i.e. sprites)
        self.sprite_list.draw()
        
        # draw HUD distance line between lander and spacecraft
        if self.reference != self.planet and self.orbit_zoom:
            x1 = self.lander._x
            y1 = self.lander._y
            x2 = self.spacecraft._x
            y2 = self.spacecraft._y
            xt = 0.5*(x1 + x2)
            yt = 0.5*(y1 + y2)
            arcade.draw_line(x1, y1, x2, y2, GREEN, 1)
            dist = math.sqrt((self.lander.x - self.spacecraft.x)**2 + (self.lander.y - self.spacecraft.y)**2)
            txt = '{:.1f}m'.format(dist)
            if dist > 1e3:
                txt = '{:.1f}km'.format(1e-3*dist)
            arcade.draw_text(txt, xt, yt, GREEN,  HUD_FONT_SIZE, font_name=HUD_FONT)
                
        arcade.draw_text("Press T to pause and view the tutorial, ESC to exit!", 10, SCREEN_HEIGHT - 30, arcade.color.WHITE, 14)
        
        # draw information
        color = arcade.color.GRAY
        size  = 12
        y_offset = SCREEN_HEIGHT - 20  # Starting position from the top of the screen
        x_offset = SCREEN_WIDTH - 220

        arcade.draw_text(f"Scale: {self.scaleFactor():.2e}", x_offset, y_offset, color, size)
        y_offset -= 20
        arcade.draw_text(f"Time Factor: {self.time_factor:.2f}", x_offset, y_offset, color, size)
        y_offset -= 20
        arcade.draw_text(f"Time: {self.sim_time:.2f}", x_offset, y_offset, color, size)
        y_offset -= 20
        arcade.draw_text("Fuel", x_offset, y_offset, color, size)
        y_offset -= 12
        self.draw_fuel_bar(x_offset, y_offset)
        y_offset -= 20
        arcade.draw_text(f"Reference: {self.reference.name}", x_offset, y_offset, color, size)
        y_offset -= 20
        arcade.draw_text(f"Control: {self.control_craft.name}", x_offset, y_offset, color, size)
        y_offset -= 20
        self.draw_control_craft_icon(x_offset + 20 + 10, y_offset)
        y_offset -= 40
        # arcade.draw_text(f"Alpha: {math.degrees(self.control_craft.alpha):.2f}°", x_offset, y_offset, color, size)
    
        if self.show_crashed:
            arcade.draw_text("Crashed!", 350, 300, arcade.color.RED, 40)

    def draw_control_craft_icon(self, x, y):
        if self.control_craft.name == "Moon":
            arcade.draw_circle_filled(x, y, 20 // 2, arcade.color.LIGHT_GRAY)  # A simple circle for the Moon
        elif self.control_craft.name == "Main Craft":
            texture = arcade.load_texture(IMAGE_PREFIX + "craft_01.png")
            scale = .06
            arcade.draw_texture_rect(
                texture,
                arcade.XYWH(x, y, texture.width, texture.height).scale(scale)
            )
        elif self.control_craft.name == "Lander":
            texture = arcade.load_texture(IMAGE_PREFIX + "lander_01.png")
            scale = .05
            arcade.draw_texture_rect(
                texture,
                arcade.XYWH(x, y, texture.width, texture.height).scale(scale)
            )
    
    def draw_fuel_bar(self, x, y):
        fuel_percentage = 100 * self.control_craft.fuel / self.control_craft.max_fuel
        bar_width = 150
        bar_height = 8
        arcade.draw_lrbt_rectangle_filled(x, x + bar_width, y, y + bar_height, arcade.color.DARK_GRAY)
        fuel_color = arcade.color.GREEN if fuel_percentage > 20 else arcade.color.RED
        arcade.draw_lrbt_rectangle_filled(x, x + (fuel_percentage / 100) * bar_width, y, y + bar_height, fuel_color)
        arcade.draw_text(f"{fuel_percentage:.1f}%", x + bar_width + 10, y - bar_height // 2, arcade.color.WHITE, HUD_FONT_SIZE)


    def on_update(self, delta_time):
        if self.paused:  
            return  
        if  self.game_running:       
            self.lander.autoPilot(delta_time * self.time_factor)
            phys_time_1 = datetime.datetime.now()
            N = max(1, int(delta_time * self.time_factor / self.sim_dt))
            dt = self.time_factor * delta_time / N                
            for i in range(N):
                for body in self.body_list:
                    # for space crafts compute the total mass
                    if isinstance(body, Spacecraft):
                        body.mass = body.dry_mass + body.fuel
                        if body.docked_to is not None:
                            body.docked_to.mass += body.mass
                for body in self.body_list:
                    body.update_physics(dt)
                for body in self.body_list:
                    body.update_metrics()
            self.sprite_list.update(delta_time)
            self.sim_time += delta_time * self.time_factor
            phys_time_2 = datetime.datetime.now()
            if (phys_time_2 - phys_time_1).total_seconds() > 0.9 * delta_time:
                self.time_factor *= 0.1     
            if (self.lander.crashed == True or self.spacecraft.crashed == True):
                self.game_running = False
                self.show_crashed = True
                self.game_over()
                
    def game_over(self):
        # print("Game Over!")
        if self.show_crashed:
            arcade.schedule_once(self.switch_to_game_over_view, 4)
        else:
            game_over_view = GameOver() 
            self.window.show_view(game_over_view) 

    def switch_to_game_over_view(self, delta_time):
        self.show_crashed = False 
        game_over_view = GameOver() 
        self.window.show_view(game_over_view) 

    def on_key_press(self, key, modifiers):
        if key == arcade.key.T:
            self.paused = True 
            tutorial = TutorialView(self)  # Pass `self` to resume later
            self.window.show_view(tutorial)
        
        elif key == arcade.key.ESCAPE:
            self.game_running = False 
            self.game_over()

        elif key == arcade.key.SPACE:
            # find index of current reference body
            i = self.body_list.index(self.reference)
            self.reference.scale_factor = self.scale_factor
            # switch to next body in list            
            i = (i + 1) % len(self.body_list)
            while self.body_list[i].docked_to is not None:
                i = (i + 1) % len(self.body_list)
            self.reference = self.body_list[i]
            self.scale_factor = self.reference.scale_factor
            
        elif key == arcade.key.NUM_ADD:
            if modifiers & arcade.key.MOD_SHIFT:
                self.time_factor *= 10
            else:
                if self.reference == self.planet:
                    self.reference.scale_factor *= 1.1
                elif self.orbit_zoom:
                    self.reference.scale_factor *= 1.1
                else:
                    self.reference.zoom_factor  *= 1.1
                
        elif key == arcade.key.NUM_SUBTRACT:
            if modifiers & arcade.key.MOD_SHIFT:
                self.time_factor /= 10
                self.time_factor = max(0.001, self.time_factor)
            else:
                if self.reference == self.planet:
                    self.reference.scale_factor /= 1.1
                elif self.orbit_zoom:
                    self.reference.scale_factor /= 1.1
                else:
                    self.reference.zoom_factor  /= 1.1
                
        elif key == arcade.key.UP:
            if self.control_craft.auto_pilot:
                self.control_craft.tgt_vvel = min(20.0, self.control_craft.tgt_vvel + 1)
            else:
                if modifiers & arcade.key.MOD_SHIFT:
                    self.control_craft.increase_thrust(10)
                else:
                    self.control_craft.increase_thrust(1)
            
        elif key == arcade.key.DOWN:
            if self.control_craft.auto_pilot:
                self.control_craft.tgt_vvel = max(-20.0, self.control_craft.tgt_vvel - 1)
            else:
                if modifiers & arcade.key.MOD_SHIFT:
                    self.control_craft.decrease_thrust(10)
                else:
                    self.control_craft.decrease_thrust(1)
                
        elif key == arcade.key.PERIOD:
            self.control_craft.increase_thrust(0.01)
            
        elif key == arcade.key.COMMA:
            self.control_craft.decrease_thrust(0.01)
            
        elif key == arcade.key.GREATER:
            self.control_craft.increase_thrust(0.1)
            
        elif key == arcade.key.LESS:
            self.control_craft.decrease_thrust(0.1)
            
        elif key == arcade.key.LEFT:
            w  = self.control_craft.omega
            w -= 0.01*2*math.pi
            w  = min(math.pi, max(-math.pi, w))
            self.control_craft.omega = w
            
        elif key == arcade.key.RIGHT:
            w  = self.control_craft.omega
            w += 0.01*2*math.pi
            w  = min(math.pi, max(-math.pi, w))
            self.control_craft.omega = w
            
        elif key == arcade.key.NUM_MULTIPLY:
            if self.control_craft.trajectory is None:
                self.control_craft.compute_orbit_2d(self.planet)
            else:
                self.control_craft.trajectory = None
                
        elif key == arcade.key.NUM_ENTER:
            self.orbit_zoom = not self.orbit_zoom
            
        elif key == arcade.key.NUM_0:
            if self.orbit_zoom:
                if self.lander.docked_to is None:
                    self.lander.dock()
                else:
                    self.lander.undock()
                
        elif key == arcade.key.NUM_1:
            self.control_craft = self.spacecraft
            
        elif key == arcade.key.NUM_2:
            self.control_craft = self.lander
            
        elif key == arcade.key.D:
            if self.control_craft == self.spacecraft:
                if self.control_craft.target is None:
                    self.control_craft.target = self.lander
                else:
                    self.control_craft.target = None
            elif self.control_craft == self.lander:
                if self.control_craft.target is None:
                    self.control_craft.target = self.spacecraft
                else:
                    self.control_craft.target = None
                    
        elif key == arcade.key.A:
            if self.control_craft == self.lander:
                self.control_craft.auto_pilot = not self.control_craft.auto_pilot
                if self.control_craft.auto_pilot:
                    self.control_craft.resetAutoPilot()
                    self.control_craft.tgt_vvel = min(20.0, max(-20.0, round(self.control_craft.verticalVelocity())))
                    
        #
        # debug final approach
        #
        elif key == arcade.key.F5 and modifiers & arcade.key.MOD_SHIFT:
            L = self.lander
            self.control_craft = L
            self.reference = L
            L.u = 0
            L.v = 0
            L.alpha = 0
            L.omega = 0
            L.acc_omega = 0
            L.acc_x = 0
            L.acc_y = 0
            L.x = 0
            L.y = -(self.planet.radius + 10e3)
            L.docked_to = None
            L.dock()


class TutorialView(arcade.View):
    def __init__(self, game_view):
        super().__init__()
        self.game_view = game_view  

    def on_draw(self):
        self.clear()
        
        arcade.draw_text("Lander Challenge", SCREEN_WIDTH / 2, SCREEN_HEIGHT - 80,
                         arcade.color.YELLOW, font_size=30, anchor_x="center")
        arcade.draw_text("enGits Lander Challenge", SCREEN_WIDTH / 2, SCREEN_HEIGHT - 120,
                         arcade.color.WHITE, font_size=20, anchor_x="center")

        tutorial_text = [
            "Press ESC to return",
            "Press ENTER to restart",
            "",
            "## Arrow Keys",
            "- LEFT: Rotate left",
            "- RIGHT: Rotate right",
            "- UP: Increase thrust by 1%",
            "- DOWN: Decrease thrust by 1%",
            "- SHIFT + UP: Increase thrust by 10%",
            "- SHIFT + DOWN: Decrease thrust by 10%",
            "",
            "## NUM-PAD",
            "- 1: Switch control to main craft",
            "- 2: Switch control to lander (if not docked)",
            "- 0: Undock/dock",
            "- ENTER: Switch to docking zoom behaviour (required for docking actions)",
            "- *: Compute/update trajectory (if Moon is reference object)",
            "- +: Zoom in",
            "- -: Zoom out",
            "- SHIFT +: Increase simulation speed",
            "- SHIFT -: Decrease simulation speed",
            "",
            "## Other Keys",
            "- SPACE: Cycle through reference objects",
            "- ,: Decrease thrust by 0.01%",
            "- <: Decrease thrust by 0.1%",
            "- .: Increase thrust by 0.01%",
            "- >: Increase thrust by 0.1%",
            "- A: Autopilot for landing below 5 km",

            "- SHIFT + F5: Put lander 10km above lunar surface (for debugging)",
            "",
            "Press ESC to return"
        ]

        y_position = SCREEN_HEIGHT - 160
        for line in tutorial_text:
            arcade.draw_text(line, SCREEN_WIDTH / 2, y_position,
                             arcade.color.WHITE if "##" not in line else arcade.color.YELLOW,
                             font_size=16 if "##" not in line else 18,
                             anchor_x="center")
            y_position -= 30 

    def on_key_press(self, key, modifiers):
        if key == arcade.key.ESCAPE:
            self.game_view.paused = False  
            self.window.show_view(self.game_view)  

        elif key == arcade.key.ENTER or key == arcade.key.NUM_ENTER:
            new_game = OrbitGame() 
            self.window.show_view(new_game)


class MainMenu(arcade.View):
    def __init__(self):
        super().__init__()
        arcade.set_background_color(arcade.color.BLACK)
        self.stars = [(random.randint(0, SCREEN_WIDTH), random.randint(0, SCREEN_HEIGHT)) for _ in range(100)]  

    def on_draw(self):
        self.clear()
        for star in self.stars:
            arcade.draw_circle_filled(star[0], star[1], 1, arcade.color.WHITE) 
        texture = arcade.load_texture(IMAGE_PREFIX + "logo.png")
        scale = .5
        arcade.draw_texture_rect(
            texture,
            arcade.XYWH(SCREEN_WIDTH / 2, SCREEN_HEIGHT / 2 + 250, texture.width, texture.height).scale(scale)
        )
        arcade.draw_text("Orbit Game", SCREEN_WIDTH / 2, SCREEN_HEIGHT / 2 + 100,
                         arcade.color.BLEU_DE_FRANCE, font_size=40, anchor_x="center")
        arcade.draw_text("Press ENTER to Start", SCREEN_WIDTH / 2, SCREEN_HEIGHT / 2 + 20,
                         arcade.color.WHITE, font_size=20, anchor_x="center")
        arcade.draw_text("Press T for Tutorial", SCREEN_WIDTH / 2, SCREEN_HEIGHT / 2 - 20,
                         arcade.color.YELLOW, font_size=20, anchor_x="center")
        arcade.draw_text("Press SHIFT-Q to Quit", SCREEN_WIDTH / 2, SCREEN_HEIGHT / 2 - 60,
                         arcade.color.WHITE, font_size=20, anchor_x="center")

    def on_key_press(self, key, modifiers):
        if key == arcade.key.ENTER or key == arcade.key.NUM_ENTER:
            self.start_game()
        elif key == arcade.key.T: 
            tutorial = TutorialView(self)
            self.window.show_view(tutorial)
        elif key == arcade.key.Q and modifiers & arcade.key.MOD_SHIFT:
            arcade.exit()

    def start_game(self):
        # print("Game Starts!")
        menu = LevelMenu()  
        self.window.show_view(menu)  


class LevelMenu(arcade.View):
    def __init__(self):
        super().__init__()
        arcade.set_background_color(arcade.color.BLACK)
        self.stars = [(random.randint(0, SCREEN_WIDTH), random.randint(0, SCREEN_HEIGHT)) for _ in range(100)]  

    def on_draw(self):
        self.clear()
        for star in self.stars:
            arcade.draw_circle_filled(star[0], star[1], 1, arcade.color.WHITE) 

        arcade.draw_text("1. Full Mission", SCREEN_WIDTH / 2, SCREEN_HEIGHT / 2 + 60,
                         arcade.color.WHITE, font_size=20, anchor_x="center")
        arcade.draw_text("2. Final Approach", SCREEN_WIDTH / 2, SCREEN_HEIGHT / 2 + 20,
                         arcade.color.WHITE, font_size=20, anchor_x="center")
        arcade.draw_text("3. Docking", SCREEN_WIDTH / 2, SCREEN_HEIGHT / 2 - 20,
                         arcade.color.WHITE, font_size=20, anchor_x="center")
        arcade.draw_text("Press SHIFT-Q to Quit", SCREEN_WIDTH / 2, SCREEN_HEIGHT / 2 - 60,
                         arcade.color.YELLOW, font_size=20, anchor_x="center")

    def on_key_press(self, key, modifiers):
        if key == arcade.key.KEY_1 or key == arcade.key.NUM_1:
            self.start_game(1)
        elif key == arcade.key.KEY_2 or key == arcade.key.NUM_2: 
            self.start_game(2)
        elif key == arcade.key.KEY_3 or key == arcade.key.NUM_3: 
            self.start_game(3)
        elif key == arcade.key.Q and modifiers & arcade.key.MOD_SHIFT:
            arcade.exit()

    def start_game(self, level):
        # print("Game Starts!")
        game = OrbitGame(level)  
        self.window.show_view(game)  


class GameOver(arcade.View):
    def __init__(self):
        super().__init__()
        arcade.set_background_color(arcade.color.BLACK)
        self.stars = [(random.randint(0, SCREEN_WIDTH), random.randint(0, SCREEN_HEIGHT)) for _ in range(100)]  # 100 stars

    def on_draw(self):
        self.clear()
        for star in self.stars:
            arcade.draw_circle_filled(star[0], star[1], 1, arcade.color.WHITE)
        arcade.draw_text("Game Over", SCREEN_WIDTH / 2, SCREEN_HEIGHT / 2 + 100,
                         arcade.color.RED, font_size=40, anchor_x="center")
        arcade.draw_text("Press ENTER to Restart", SCREEN_WIDTH / 2, SCREEN_HEIGHT / 2 + 20,
                         arcade.color.WHITE, font_size=20, anchor_x="center")
        arcade.draw_text("Press T for Tutorial", SCREEN_WIDTH / 2, SCREEN_HEIGHT / 2 - 20,
                         arcade.color.YELLOW, font_size=20, anchor_x="center")
        arcade.draw_text("Press SHIFT-Q", SCREEN_WIDTH / 2, SCREEN_HEIGHT / 2 - 60,
                         arcade.color.WHITE, font_size=20, anchor_x="center")

    def on_key_press(self, key, modifiers):
        if key == arcade.key.ENTER or key == arcade.key.NUM_ENTER:
            self.restart_game() 
        elif key == arcade.key.T: 
            tutorial = TutorialView(self)
            self.window.show_view(tutorial)
        elif key == arcade.key.Q and modifiers & arcade.key.MOD_SHIFT:
            arcade.exit()

    def restart_game(self):
        # print("Game Restarts!")
        self.window.clear()
        game = MainMenu() 
        self.window.show_view(game) 


def main():
    window = arcade.Window(SCREEN_WIDTH, SCREEN_HEIGHT, SCREEN_TITLE)
    menu = MainMenu()  
    window.show_view(menu) 
    arcade.run() 


if __name__ == "__main__":
    main()