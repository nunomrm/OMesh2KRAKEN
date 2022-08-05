      program ramctr
c
c     Version 1.4
c
c     This code creates contour TL plots from the output of RAM. 
c     The Postscript plot file is ramctr.ps. The subroutines emulate 
c     subroutines from the DISSPLA graphics package.
c
      real tl(1000,1000),xb(4),yb(4),rb(1000),zb(1000),
     >   xc(4),yc(4)
      open(unit=1,status='old',file='ram.in')
      open(unit=2,status='old',file='tl.grid',form='unformatted')
c
      read(1,*)
      read(1,*)
      read(1,*)xmax
      read(1,*)zdum,dz,ndz,ymax
      read(1,*)
      xmax=xmax*0.001
c
      nb=1
    1 read(1,*)rb(nb),zb(nb)
      rb(nb)=rb(nb)*0.001
      if(rb(nb).lt.0.0)go to 2
      nb=nb+1
      go to 1
    2 rb(nb)=xmax
      zb(nb)=zb(nb-1)
c
      write(*,*)' '
      write(*,*)'   enter delr (km) and delz (m)'
      write(*,*)' '
      read(*,*)delx,dely
c
      write(*,*)' '
      write(*,*)'   enter tlmin, tlmax, deltl (dB)'
      write(*,*)' '
      read(*,*)tlmin,tlmax,deltl
      write(*,*)' '
c
      read(2)ny
      nx=1
    3 read(2,end=4)(tl(j,nx),j=1,ny)
      nx=nx+1
      go to 3
c
    4 nx=nx-1
c
      call comprs
      call area2d(6.0,4.0)
      call xname('Range (km)',10)
      call yname('Depth (m)',9)
      call frame
      call graf(0.0,delx,xmax,ymax,-dely,0.0)
c
      dx=xmax/float(nx)
      dy=ymax/float(ny)
c
      tlc=tlmin
    5 call contr(tl,tlc,xc,yc,nc,nx,ny,xmin,ymin,dx,dy)
      tlc=tlc+deltl
      if(tlc.lt.tlmax+0.5*deltl)go to 5
c
      call thkcrv(0.02)
      call curve(rb,zb,nb,0)
c
      call endpl(0)
      call donepl
c
      stop
      end
c
      subroutine comprs
      common w,h,xx0,yy0,xx1,yy1,ax,bx,ay,by,ifrm,idsh
c
      open(unit=40,status='unknown',file='ramctr.ps')
      write(40,'(a)')'%!'
      write(40,*)'/center {'
      write(40,*)'    dup stringwidth pop'
      write(40,*)'    2 div neg 0 rmoveto'
      write(40,*)'} bind def'
      ifrm=0
      idsh=0
c
      return
      end
c
      subroutine area2d(ww,hh)
      common w,h,xx0,yy0,xx1,yy1,ax,bx,ay,by,ifrm,idsh
c
      w=ww
      h=hh
      xx0=(612.0-72.0*w+36.0)/2.0
      yy0=(792.0-72.0*h+36.0)/2.0
      xx1=xx0+72.0*w
      yy1=yy0+72.0*h
c
      return
      end
c
      subroutine frame
      common w,h,xx0,yy0,xx1,yy1,ax,bx,ay,by,ifrm,idsh
      ifrm=1
      return
      end
c
      subroutine graf(xmin,delx,xmax,ymin,dely,ymax)
      common w,h,xx0,yy0,xx1,yy1,ax,bx,ay,by,ifrm,idsh
      eps=0.000001
c
      write(40,*)'/Times-Roman findfont'
      write(40,*)'14 scalefont'
      write(40,*)'setfont'
c
      ax=xx0-w*72.0*xmin/(xmax-xmin)
      bx=w*72.0/(xmax-xmin)
      ay=yy0-h*72.0*ymin/(ymax-ymin)
      by=h*72.0/(ymax-ymin)
c
      mx=ifix(0.5+(xmax-xmin)/delx)
      xtick=xmin
      do 1 i=1,mx+1
      xx=ax+bx*xtick
      yy=yy0
      write(40,*)xx,yy,' moveto'
      yy=yy0-9.0
      write(40,*)xx,yy,' lineto'
      write(40,*)'stroke'
      yy=yy0-22.0
      write(40,*)xx,yy,' moveto'
      jtick=ifix(alog10(abs(xtick)+eps)+eps)
      jtick=max(0,jtick)
      if(xtick.gt.-eps)then
      itick=ifix(xtick+0.5)
      else
      itick=ifix(xtick-0.5)
      jtick=jtick+1
      end if
      if(jtick.eq.0)write(40,3)itick
      if(jtick.eq.1)write(40,4)itick
      if(jtick.eq.2)write(40,5)itick
      if(jtick.eq.3)write(40,6)itick
      if(jtick.eq.4)write(40,7)itick
      xtick=xtick+delx
    1 continue
c
      my=ifix(0.5+(ymax-ymin)/dely)
      ytick=ymin
      do 2 i=1,my+1
      xx=xx0
      yy=ay+by*ytick
      write(40,*)xx,yy,' moveto'
      xx=xx0-9.0
      write(40,*)xx,yy,' lineto'
      write(40,*)'stroke'
      xx=xx0-12.0
      write(40,*)xx,yy,' moveto'
      jtick=ifix(alog10(abs(ytick)+eps)+eps)
      jtick=max(0,jtick)
      if(ytick.gt.-eps)then
      itick=ifix(ytick+0.5)
      else
      itick=ifix(ytick-0.5)
      jtick=jtick+1
      end if
      write(40,*)'90 rotate'
      if(jtick.eq.0)write(40,3)itick
      if(jtick.eq.1)write(40,4)itick
      if(jtick.eq.2)write(40,5)itick
      if(jtick.eq.3)write(40,6)itick
      if(jtick.eq.4)write(40,7)itick
      write(40,*)'270 rotate'
      ytick=ytick+dely
    2 continue
c
    3 format('(',i1,')',' center show')
    4 format('(',i2,')',' center show')
    5 format('(',i3,')',' center show')
    6 format('(',i4,')',' center show')
    7 format('(',i5,')',' center show')
c
      return
      end
c
      subroutine curve(x,y,n,icrv)
      real x(5000),y(5000)
      common w,h,xx0,yy0,xx1,yy1,ax,bx,ay,by,ifrm,idsh
c
c     Solid curve.
c
      if(idsh.eq.0)then
      xx=ax+bx*x(1)
      yy=ay+by*y(1)
      write(40,*)xx,yy,' moveto'
      do 1 i=2,n
      xx=ax+bx*x(i)
      yy=ay+by*y(i)
      write(40,*)xx,yy,' lineto'
    1 continue
      write(40,*)'stroke'
      return
      end if
c
c     Dashed curve.
c
      s=0.0
      dels=10.0
      iblck=1
      xx=ax+bx*x(1)
      yy=ay+by*y(1)
      write(40,*)xx,yy,' moveto'
      xxold=xx
      yyold=yy
c
      do 3 i=2,n
      xx=ax+bx*x(i)
      yy=ay+by*y(i)
c
    2 dx=xx-xxold
      dy=yy-yyold
      ds=sqrt(dx**2+dy**2)
      s=s+ds
c
      if(s.le.dels)then
      xxold=xx
      yyold=yy
      if(iblck.eq.1)write(40,*)xx,yy,' lineto'
      go to 3
      end if
c
      frac=1.0-(s-dels)/ds
      xxold=xxold+frac*dx
      yyold=yyold+frac*dy
c
      if(iblck.eq.1)then
      write(40,*)xxold,yyold,' lineto'
      write(40,*)'stroke'
      iblck=0
      dels=3.0
c
      else
      iblck=1
      dels=10.0
      end if
c
      write(40,*)xxold,yyold,' moveto'
      s=0.0
      go to 2
c
    3 continue
c
      if(iblck.eq.1)write(40,*)'stroke'
c
      return
      end
c
      subroutine xname(xchar,n)
      character*50 xout,xchar
      common w,h,xx0,yy0,xx1,yy1,ax,bx,ay,by,ifrm,idsh
c
      xx=0.5*(xx0+xx1)
      yy=yy0-44.0
      xout='('//xchar(1:n)//') center show'
      write(40,*)'/Times-Roman findfont'
      write(40,*)'18 scalefont'
      write(40,*)'setfont'
      write(40,*)xx, yy,' moveto'
      write(40,*)xout
c
      return
      end
c
      subroutine yname(ychar,n)
      character*50 yout,ychar
      common w,h,xx0,yy0,xx1,yy1,ax,bx,ay,by,ifrm,idsh
c
      xx=xx0-32.0
      yy=0.5*(yy0+yy1)
      yout='('//ychar(1:n)//') center show'
      write(40,*)'/Times-Roman findfont'
      write(40,*)'18 scalefont'
      write(40,*)'setfont'
      write(40,*)xx, yy,' moveto'
      write(40,*)'90 rotate'
      write(40,*)yout
      write(40,*)'270 rotate'
c
      return
      end
c
      subroutine dash
      common w,h,xx0,yy0,xx1,yy1,ax,bx,ay,by,ifrm,idsh
      idsh=idsh-1
      return
      end
c
      subroutine shade(xb,yb,idum1,dum1,dum2,idum2,idum3,idum4)
      real xb(4),yb(4)
      common w,h,xx0,yy0,xx1,yy1,ax,bx,ay,by,ifrm,idsh
c
      xx=ax+bx*xb(1)
      yy=ay+by*yb(1)
      write(40,*)xx,yy,' moveto'
      xx=ax+bx*xb(2)
      yy=ay+by*yb(2)
      write(40,*)xx,yy,' lineto'
      xx=ax+bx*xb(3)
      yy=ay+by*yb(3)
      write(40,*)xx,yy,' lineto'
      xx=ax+bx*xb(4)
      yy=ay+by*yb(4)
      write(40,*)xx,yy,' lineto'
      write(40,*)'closepath'
      write(40,*)'fill'
c
      return
      end
c
      subroutine thkcrv(thknss)
      write(40,*)thknss*72.0,' setlinewidth'
      return
      end
c
      subroutine hwhsi(hue,sat,brt)
      write(40,*)hue,sat,brt,' sethsbcolor'
      return
      end
c
c     Construct contour curves.
c
      subroutine contr(f,f0,xc,yc,nc,nx,ny,xmin,ymin,dx,dy)
      common w,h,xx0,yy0,xx1,yy1,ax,bx,ay,by,ifrm,idsh
      real f(1000,1000),xc(4),yc(4)
c
c     Draw the contour curve segments inside each rectangle.
c
      do 2 i=1,ny-1
      yy=float(i)*dy
c
      do 1 j=1,nx-1
      xx=float(j)*dx
      nc=0
c
      f1=f(i,j)
      f2=f(i+1,j)
      f3=f(i,j+1)
      f4=f(i+1,j+1)
c
      det12=(f1-f0)*(f2-f0)
      det24=(f2-f0)*(f4-f0)
      det13=(f1-f0)*(f3-f0)
      det34=(f3-f0)*(f4-f0)
c
      if(det13.le.0.0)then
      nc=nc+1
      xc(nc)=xx+dx*(f0-f1)/(f3-f1)
      yc(nc)=yy
      end if
c
      if(det34.le.0.0)then
      nc=nc+1
      xc(nc)=xx+dx
      yc(nc)=yy+dy*(f0-f3)/(f4-f3)
      end if
c
      if(det24.le.0.0)then
      nc=nc+1
      xc(nc)=xx+dx*(f0-f2)/(f4-f2)
      yc(nc)=yy+dy
      end if
c
      if(det12.le.0.0)then
      nc=nc+1
      xc(nc)=xx
      yc(nc)=yy+dy*(f0-f1)/(f2-f1)
      end if
c
      if(nc.ge.2)call cnect(f,f1,f2,f3,f4,xc,yc,nc,dx,dy,xx,yy)
c
    1 continue
    2 continue
c
      return
      end
c
c     Connect the contour points.
c
      subroutine cnect(f,f1,f2,f3,f4,xc,yc,nc,dx,dy,xx,yy)
      common w,h,xx0,yy0,xx1,yy1,ax,bx,ay,by,ifrm,idsh
      real f(1000,1000),xc(4),yc(4)
c
      if(nc.le.3)then
      xxc=ax+bx*xc(1)
      yyc=ay+by*yc(1)
      write(40,*)xxc,yyc,' moveto'
      xxc=ax+bx*xc(2)
      yyc=ay+by*yc(2)
      write(40,*)xxc,yyc,' lineto'
      write(40,*)' stroke'
      return
      end if
c
c     The saddle point case. 
c     Approximate f(x,y) by a*x + b*y + c*x*y + d. 
c     The mappings x = x(f0,y) and y = y(f0,x) are
c     singular at y = -a/c and x = -b/c. 
c
      b=(xx+dx)*(f2-f1)-xx*(f4-f3)
      c=f1-f2-f3+f4
      xsng=-b/c
      det=(xc(1)-xsng)*(xc(2)-xsng)
c
      if(det.lt.0.0)then
      xxc=ax+bx*xc(1)
      yyc=ay+by*yc(1)
      write(40,*)xxc,yyc,' moveto'
      xxc=ax+bx*xc(4)
      yyc=ay+by*yc(4)
      write(40,*)xxc,yyc,' lineto'
      write(40,*)' stroke'
      xxc=ax+bx*xc(2)
      yyc=ay+by*yc(2)
      write(40,*)xxc,yyc,' moveto'
      xxc=ax+bx*xc(3)
      yyc=ay+by*yc(3)
      write(40,*)xxc,yyc,' lineto'
      write(40,*)' stroke'
c
      else
      xxc=ax+bx*xc(1)
      yyc=ay+by*yc(1)
      write(40,*)xxc,yyc,' moveto'
      xxc=ax+bx*xc(2)
      yyc=ay+by*yc(2)
      write(40,*)xxc,yyc,' lineto'
      write(40,*)' stroke'
      xxc=ax+bx*xc(3)
      yyc=ay+by*yc(3)
      write(40,*)xxc,yyc,' moveto'
      xxc=ax+bx*xc(4)
      yyc=ay+by*yc(4)
      write(40,*)xxc,yyc,' lineto'
      write(40,*)' stroke'
      end if
c
      return
      end
c
      subroutine endpl(iplot)
      common w,h,xx0,yy0,xx1,yy1,ax,bx,ay,by,ifrm,idsh
c
      if(ifrm.eq.1)then
      write(40,*)'2.5 setlinewidth'
      write(40,*)xx0,yy0,' moveto'
      write(40,*)xx1,yy0,' lineto'
      write(40,*)xx1,yy1,' lineto'
      write(40,*)xx0,yy1,' lineto'
      write(40,*)xx0,yy0,' lineto'
      write(40,*)xx1,yy0,' lineto'
      write(40,*)'stroke'
      write(40,*)'1 setlinewidth'
      end if
c
      return
      end
c
      subroutine donepl
      write(40,*)'showpage'
      return
      end
