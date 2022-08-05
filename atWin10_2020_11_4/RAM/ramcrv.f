      program ramcrv
c
c     Version 1.4
c
c     This code creates TL curves from the output of RAM. 
c     The Postscript plot file is ramclr.ps. The subroutines 
c     emulate subroutines from the DISSPLA graphics package. 
c     The output of RAM is displayed as a solid curve. A 
c     dashed curve is included for comparison if the file 
c     tl.dash is present. This file is in the same format as 
c     the output file tl.line of RAM.  
c
      real x(5000),y(5000)
      open(unit=1,status='old',file='tl.line')
c
      write(*,*)' '
      write(*,*)'   enter rmin, rmax, delr (km)'
      write(*,*)' '
      read(*,*)xmin,xmax,delx
c
      write(*,*)' '
      write(*,*)'   enter tlmin, tlmax, deltl'
      write(*,*)' '
      read(*,*)ymax,ymin,dely
      write(*,*)' '
c
      call comprs
      call area2d(6.0,4.0)
      call xname('Range (km)',10)
      call yname('Loss (dB re 1m)',15)
      call frame
      call graf(xmin,delx,xmax,ymin,-dely,ymax)
c
      n=0
    1 read(1,*,end=2)r,tl
      r=r*0.001
      if(r.lt.xmin)go to 1
      if(r.gt.xmax)go to 2
      n=n+1
      x(n)=r
      tl=amin1(tl,ymin)
      tl=amax1(tl,ymax)
      y(n)=tl
      go to 1
c
    2 call curve(x,y,n,0)
c
      open(unit=2,status='old',file='tl.dash',err=5)
c
      n=0
    3 read(2,*,end=4)r,tl
      r=r*0.001
      if(r.lt.xmin)go to 3
      if(r.gt.xmax)go to 4
      n=n+1
      x(n)=r
      tl=amin1(tl,ymin)
      tl=amax1(tl,ymax)
      y(n)=tl
      go to 3
c
    4 call dash
      call curve(x,y,n,0)
c
    5 call endpl(0)
      call donepl
c
      stop
      end
c
      subroutine comprs
      common w,h,xx0,yy0,xx1,yy1,ax,bx,ay,by,ifrm,idsh
c
      open(unit=40,status='unknown',file='ramcrv.ps')
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
