import arcade
import math
import sys
import datetime
import numpy as np


SCREEN_WIDTH   = 1200
SCREEN_HEIGHT  = 1200
SCREEN_TITLE   = 'enGits Lunar Lander Challenge'
NUM_MOUNTAINS  = 1000
ALT0           = 2.0
FEAT_THRESHOLD = 1000
MOON_COLOR     = arcade.color.DARK_GRAY


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
            # if _x < -FEAT_THRESHOLD or _x > SCREEN_WIDTH + FEAT_THRESHOLD:
            #     return
            # if _y < -FEAT_THRESHOLD or _y > SCREEN_HEIGHT + FEAT_THRESHOLD:
            #     return
            pts.append((_x, _y))
        arcade.draw_polygon_filled(pts, self.color)

    
class Mountain(SurfaceFeature):
    
    def __init__(self, planet, height, width, middle_width=None, theta=0):  
        super().__init__(planet, theta)
        self.color  = MOON_COLOR
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
        if self.sprite is None:
            self.draw_high_res_circle(self._x, self._y, self.game.scaleFactor()*self.radius, self.color, num_segments=150)
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
            self.sprite.append_texture(arcade.load_texture(_file))
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
        self.crashed      = False

    def update_metrics(self):
        r   = math.sqrt(self.x**2 + self.y**2)
        R   = self.game.planet.radius
        alt = r - R
        if alt < 10e3 and self.game.reference == self:
            self.game.scale_factor = min(0.25*min(SCREEN_WIDTH, SCREEN_HEIGHT) / alt, self.game.real_zoom_factor)
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
        grn = (0, 200, 0)
        #
        # draw HUD if required
        #
        if self == self.game.control_craft: # and self == self.game.reference:
            fs  = 10
            vel = math.sqrt(self.u**2 + self.v**2)
            alt = math.sqrt(self.x**2 + self.y**2) - self.game.planet.radius - ALT0
            x   = self._x  + 20
            y   = self._y
            fn  = 'Courier New'
            col = grn
            #
            if alt > 5.0e3:
                if self == self.game.lander:
                    col = arcade.color.ORANGE
                arcade.draw_text('km/h {:.2f}'.format(vel*3.6), x, y + fs + 2, col, fs, font_name=fn)
                arcade.draw_text('alt  {:.2f}'.format(alt/1000), x, y, col, fs, font_name=fn)
                arcade.draw_text('fuel {:.1f}%'.format(100*self.fuel/self.max_fuel), x, y - fs - 2, col, fs, font_name=fn)
                arcade.draw_text('th.  {:.2f}%'.format(self.thrust_level), x, y - 2*fs - 4, col, fs, font_name=fn)
            else:
                if self == self.game.lander:
                    col = arcade.color.RED
                if vel < 5.0:
                    col = grn
                arcade.draw_text('m/s {:.2f}'.format(vel), x, y + fs + 2, col, fs, font_name=fn)
                arcade.draw_text('m   {:.2f}'.format(alt), x, y, col, fs, font_name=fn)
                arcade.draw_text('f.  {:.1f}%'.format(100*self.fuel/self.max_fuel), x, y - fs - 2, col, fs, font_name=fn)
                arcade.draw_text('th. {:.2f}%'.format(self.thrust_level), x, y - 2*fs - 4, col, fs, font_name=fn)
            #
            u = self.u
            v = self.v
            u_crit = 1.0
            if self.target is not None:
                u -= self.target.u
                v -= self.target.v                
                u_crit = 0.01
            u_abs = math.sqrt(u**2 + v**2)
            if u_abs > u_crit:
                dx = u / u_abs
                dy = v / u_abs
                a  = Body.safe_angle(dx, dy)
                a += self.game.reference.alpha
                dx = 100*math.cos(a)
                dy = 100*math.sin(a)
                arcade.draw_line(self._x, self._y, self._x + dx, self._y + dy, col, 1)
            
        #
        # draw orbital trajectory
        #
        if self.trajectory is not None and not self.game.real_zoom and self.game.reference == self.game.planet:
            pts = []
            for x, y in self.trajectory:
                _x = (x - self.game.reference.x) * self.game.scaleFactor() + SCREEN_WIDTH/2
                _y = (y - self.game.reference.y) * self.game.scaleFactor() + SCREEN_HEIGHT/2
                pts.append((_x, _y))
            arcade.draw_line_strip(pts, grn, 1)
            n = int(len(pts)/2)
            arcade.draw_line(pts[0][0], pts[0][1], pts[n][0], pts[n][1], grn, 1)            
            km1 = (pts[0][0] - self.game.planet._x)**2 + (pts[0][1] - self.game.planet._y)**2
            km1 = 1e-3 * (math.sqrt(km1) / self.game.scaleFactor() - self.game.planet.radius)
            km2 = (pts[n][0] - self.game.planet._x)**2 + (pts[n][1] - self.game.planet._y)**2
            km2 = 1e-3 * (math.sqrt(km2) / self.game.scaleFactor() - self.game.planet.radius)
            arcade.draw_text('{:.1f}km'.format(km1), pts[0][0], pts[0][1], grn, 10)
            arcade.draw_text('{:.1f}km'.format(km2), pts[n][0], pts[n][1], grn, 10)
    
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
            u_abs = math.sqrt(self.u**2 + self.v**2)
            if u_abs > 5.0:
                # CRASH
                # arcade.exit()
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
            if abs(alpha - beta) > 20:
                # CRASH
                # arcade.exit()
                self.crashed = True
            self.x = self.game.planet.x + dx*(self.game.planet.radius + ALT0)
            self.y = self.game.planet.y + dy*(self.game.planet.radius + ALT0)
            self.trajectory = None
            self.thrust_level = 0
            self.sprite.set_texture(0)
        if update_trajectory and self.trajectory is not None:
            self.compute_orbit_2d(self.game.planet)
        super().update_physics(dt, Fx, Fy)

    def increase_thrust(self, amount):
        self.thrust_level = min(100, self.thrust_level + amount)
        idx = int(round(4*self.thrust_level/100))
        self.sprite.set_texture(idx)
        
    def decrease_thrust(self, amount):
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
    def __init__(self):
        # super().__init__(SCREEN_WIDTH, SCREEN_HEIGHT, SCREEN_TITLE)
        super().__init__()

        self.time_factor      = 1.0
        self.sim_dt           = 0.1
        self.sim_time         = 0.0
        self.real_zoom        = False
        self.real_zoom_factor = 5.74

        arcade.set_background_color(arcade.color.BLACK)
        self.body_list = []
        self.sprite_list = arcade.SpriteList()

        self.game_running = True 

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
        self.spacecraft.y            = self.planet.radius + 110e3
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
        
    def scaleFactor(self):
        if self.real_zoom:
            return self.real_zoom_factor
        return self.scale_factor
        
    def on_draw(self):
        self.clear()
        
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
        
        # draw information
        color = arcade.color.GRAY
        size  = 12
        arcade.draw_text(f"Reference: {self.reference.name}", 10, 10, color, size)
        arcade.draw_text("Scale: {:.2e}".format(self.scaleFactor()), 10, 30, color, size)
        arcade.draw_text(f"Time factor: {self.time_factor:.2f}", 10, 50, color, size)
        arcade.draw_text(f"Time: {self.sim_time:.2f}", 10, 70, color, size)
        arcade.draw_text(f"Fuel: {100*self.control_craft.fuel/self.control_craft.max_fuel:.1f} %", 10, 90, color, size)
        arcade.draw_text(f"Control: {self.control_craft.name}", 10, 110, color, size)
        arcade.draw_text('alpha: {:.2f}deg'.format(math.degrees(self.control_craft.alpha)), 10, 130, color, size)

    def on_update(self, delta_time):
        if  self.game_running:
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
                    body.update_metrics()
            self.sprite_list.update(delta_time)
            self.sim_time += delta_time * self.time_factor
            phys_time_2 = datetime.datetime.now()
            if (phys_time_2 - phys_time_1).total_seconds() > 0.9 * delta_time:
                self.time_factor *= 0.1        
            if (self.lander.crashed == True or self.spacecraft.crashed == True):
                self.game_running = False
                self.game_over()

    def game_over(self):
        """Switch to Game Over view."""
        print("Game Over!")
        game_over_view = GameOver()
        self.window.show_view(game_over_view)

    def on_key_press(self, key, modifiers):
        if key == arcade.key.ESCAPE:
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
                self.scale_factor *= 1.1
                
        elif key == arcade.key.NUM_SUBTRACT:
            if modifiers & arcade.key.MOD_SHIFT:
                self.time_factor /= 10
                self.time_factor = max(0.001, self.time_factor)
            else:
                self.scale_factor /= 1.1
                
        elif key == arcade.key.UP:
            if modifiers & arcade.key.MOD_SHIFT:
                self.control_craft.increase_thrust(10)
            else:
                self.control_craft.increase_thrust(1)
            
        elif key == arcade.key.DOWN:
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
            self.real_zoom = not self.real_zoom
            
        elif key == arcade.key.NUM_0:
            if self.real_zoom:
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
            L.y = self.planet.radius + 10e3
            L.docked_to = None
            L.dock()


# def main():
#     game = OrbitGame()
#     arcade.run()

# if __name__ == "__main__":
#     main()

class MainMenu(arcade.View):
    def on_show(self):
        arcade.set_background_color(arcade.color.DARK_BLUE)

    def on_draw(self):
        self.clear()
        arcade.draw_text("Orbit Game", SCREEN_WIDTH / 2, SCREEN_HEIGHT / 2 + 50,
                         arcade.color.WHITE, font_size=40, anchor_x="center")
        arcade.draw_text("Press ENTER to Start", SCREEN_WIDTH / 2, SCREEN_HEIGHT / 2 - 20,
                         arcade.color.WHITE, font_size=20, anchor_x="center")

    def on_key_press(self, key, modifiers):
        if key == arcade.key.ENTER:
            self.start_game()

    def start_game(self):
        """Start the game view."""
        print("Game Starts!")
        game = OrbitGame()  
        self.window.show_view(game) 


class GameOver(arcade.View):
    def on_show(self):
        arcade.set_background_color(arcade.color.DARK_RED)

    def on_draw(self):
        self.clear()
        arcade.draw_text("Game Over", SCREEN_WIDTH / 2, SCREEN_HEIGHT / 2 + 50,
                         arcade.color.WHITE, font_size=40, anchor_x="center")
        arcade.draw_text("Press ENTER to Play Again", SCREEN_WIDTH / 2, SCREEN_HEIGHT / 2 - 20,
                         arcade.color.WHITE, font_size=20, anchor_x="center")
        arcade.draw_text("Press ESC to Quit", SCREEN_WIDTH / 2, SCREEN_HEIGHT / 2 - 60,
                         arcade.color.WHITE, font_size=20, anchor_x="center")

    def on_key_press(self, key, modifiers):
        if key == arcade.key.ENTER:
            self.restart_game()  # Restart the game
        elif key == arcade.key.ESCAPE:
            arcade.exit()  # Exit the game

    def restart_game(self):
        print("Game Restarts!")
        self.window.clear()
        game = OrbitGame() 
        self.window.show_view(game) 


def main():
    window = arcade.Window(SCREEN_WIDTH, SCREEN_HEIGHT, SCREEN_TITLE)
    menu = MainMenu()  # Create the MainMenu view
    window.show_view(menu)  # Show the menu first
    arcade.run()  # Start the event loop


if __name__ == "__main__":
    main()
