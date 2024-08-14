import pygame
import math

def redraw(screen, time, dt, score, precision, s, v, taus, phis):

  (width, height) = screen.get_size()
  scale = min(width, height)

  hell = (240,240,240)
  dunkel = (64,64,64)
  grau = (200,200,200)
  blau = (100,165,220)
  rot = (240,85,70)
  gelb = (215,210,195)
    
  screen.fill(hell)

  # define a few local functions

  def display(string, pos, centered=True):
    font = pygame.font.SysFont(None, math.ceil(scale/32))
    text = font.render(string, True, (0,0,0))
    if centered:
      w, h = font.size(string)
      screen.blit(text, (pos[0]-w//2, pos[1]-h//2))
    else:
      screen.blit(text, (pos[0], pos[1]))
      screen.blit(text, (pos[0], pos[1]))

  def plot(values, origin, color, scaling, eps=None):
    N = 850
    w = width*(0.5-0.02)
    h = height*(0.5-0.02)
    
    pygame.draw.line(screen, dunkel, origin + (0,h/2), origin - (0,h/2), 2)
    pygame.draw.line(screen, dunkel, origin, origin + (w,0), 2)
    
    if len(values) >= 2:
      xx = list(origin.x + (N-i)/N*w for i in range(N))
      yy = list(origin.y - (v*scaling/2)*h for v in values)
      pts = list(zip(xx, yy))
      pygame.draw.lines(screen, color, False, pts, width=3)

  def tachometer(value):
    radius = scale*0.2
    center = pygame.Vector2(width/4, height*0.85)
    max_val = 5

    # outer ring
    pygame.draw.circle(screen, "black", center, radius, math.ceil(radius*0.02), draw_top_left=True, draw_top_right=True)

    # left (-1)
    display(str(-max_val), center+radius*pygame.Vector2(-0.6,0))
    pygame.draw.polygon(screen, "black",
      [center + radius*pygame.Vector2(-0.75,0.01),
       center + radius*pygame.Vector2(-0.99,0.01),
       center + radius*pygame.Vector2(-0.99,-0.01),
       center + radius*pygame.Vector2(-0.75,-0.01)])
    # top (0)
    display("0", center+radius*pygame.Vector2(0,-0.6))
    pygame.draw.polygon(screen, dunkel,
      [center + radius*pygame.Vector2(0.01,-0.75),
       center + radius*pygame.Vector2(0.01,-0.99),
       center + radius*pygame.Vector2(-0.01,-0.99),
       center + radius*pygame.Vector2(-0.01,-0.75)])
    # right (1)
    display(str(max_val), center+radius*pygame.Vector2(0.6,0))
    pygame.draw.polygon(screen, "black",
      [center + radius*pygame.Vector2(0.75,0.01),
       center + radius*pygame.Vector2(0.99,0.01),
       center + radius*pygame.Vector2(0.99,-0.01),
       center + radius*pygame.Vector2(0.75,-0.01)])
    # in-between ticks
    for i in [0.25, 0.5, 0.75, 1.25, 1.5, 1.75]:
      direction = radius * pygame.Vector2(math.cos(i*math.pi/2), -math.sin(i*math.pi/2))
      orthogonal = radius * pygame.Vector2(-math.sin(i*math.pi/2), -math.cos(i*math.pi/2))
      pygame.draw.polygon(screen, dunkel,
        [center + 0.8*direction - 0.01*orthogonal,
         center + 0.9*direction - 0.01*orthogonal,
         center + 0.9*direction + 0.01*orthogonal,
         center + 0.8*direction + 0.01*orthogonal])

    # needle
    angle = (1-value/max_val)*math.pi/2
    direction = radius * pygame.Vector2(math.cos(angle), -math.sin(angle))
    orthogonal = radius * pygame.Vector2(-math.sin(angle), -math.cos(angle))
    pygame.draw.polygon(screen, rot,
      [center - 0.03*orthogonal,
       center + 0.85*direction - 0.01*orthogonal,
       center + 0.85*direction + 0.01*orthogonal,
       center + 0.03*orthogonal])
    # inner circle
    pygame.draw.circle(screen, dunkel, center, radius*0.08)

    # value
    display("{: .3f}".format(value), center + radius*pygame.Vector2(0,0.2))

  def car(s, phi):
    ground = pygame.Vector2(width/4, height*0.5)

    # moving ground
    thickness = math.ceil(scale*3e-3)
    ground_w = width*0.4
    pygame.draw.line(screen, "black", ground + (-ground_w/2,0), ground + (ground_w/2,0), thickness)
    N = 6
    for i in range(1,N+2):
      t = s % 1
      pt = ground + (ground_w*(-0.5 + (i-t)/N), 0)
      pygame.draw.line(screen, "black", pt, pt + (-scale*0.03, scale*0.03), thickness)
    # hide border of ground
    rect = pygame.Rect(ground.x-ground_w/2-scale*0.04, ground.y-1, scale*0.04, scale*0.04)
    pygame.draw.rect(screen, hell, rect)
    rect = pygame.Rect(ground.x+ground_w/2, ground.y-1, ground_w/N+2, scale*0.04)
    pygame.draw.rect(screen, hell, rect)

    # wheels
    wheel_radius = scale*0.03
    pygame.draw.circle(screen, "black", ground + (-wheel_radius*2,-wheel_radius), wheel_radius, math.ceil(wheel_radius/2))
    pygame.draw.circle(screen, "black", ground + (wheel_radius*2,-wheel_radius), wheel_radius, math.ceil(wheel_radius/2))

    # top center of car
    car_pos = ground + (0, -5*wheel_radius)

    # chassis
    w, h = 8*wheel_radius, 3.7*wheel_radius
    rect = pygame.Rect(car_pos.x-w/2, car_pos.y, w, h)
    pygame.draw.rect(screen, grau, rect, border_radius=math.ceil(scale/256))
    pygame.draw.rect(screen, "black", rect, border_radius=math.ceil(scale/256),width=thickness)

    # pendulum
    pendulum_len = scale*0.15
    pendulum_radius = scale*0.025
    pendulum_center = car_pos + pendulum_len * pygame.Vector2(math.cos(phi), -math.sin(phi))
    pygame.draw.line(screen, dunkel, car_pos, pendulum_center, 2*thickness)
    pygame.draw.circle(screen, dunkel, car_pos, 3*thickness)
    pygame.draw.circle(screen, blau, pendulum_center, pendulum_radius)
    pygame.draw.circle(screen, "black", pendulum_center, pendulum_radius, thickness)
    # light speck
    pygame.draw.circle(screen, hell, pendulum_center - (pendulum_radius/4, pendulum_radius/4), pendulum_radius*0.2)

  # do the actual drawing
  car(s, math.pi/2 + phis[-1])
  plot(phis, pygame.Vector2(width*(0.5+0.01),height/4), blau, 0.1, precision)
  tachometer(v)
  plot(taus, pygame.Vector2(width*(0.5+0.01),height*3/4), rot, 0.1)

  # FPS
  #display("FPS: {: .1f}".format(1/dt), (16,16), False)
  display("Angle: φ(t)", (3*width/4, height*0.05))
  display("Control: a(t)", (3*width/4, height*0.55))
  display("Time: t = {: .1f}".format(time), (width/4, height*0.05))


def game_over(screen, score):

  (width, height) = screen.get_size()
  scale = min(width, height)

  hell = (240,240,240)
  dunkel = (64,64,64)
  grau = (200,200,200)
  blau = (100,165,220)
  rot = (240,85,70)
  gelb = (215,210,195)

  def display(string, pos, centered=True):
    font = pygame.font.SysFont(None, math.ceil(scale/8))
    text = font.render(string, True, (0,0,0))
    if centered:
      w, h = font.size(string)
      return screen.blit(text, (pos[0]-w//2, pos[1]-h//2))
    else:
      return screen.blit(text, (pos[0], pos[1]))

  rect = display("SCORE: {:.0f}".format(score), (width/2, height/2))
  rect.scale_by_ip(1.25)
  pygame.draw.rect(screen, gelb, rect)
  pygame.draw.rect(screen, dunkel, rect, 1)
  display("SCORE: {:.0f}".format(score), (width/2, height/2))
