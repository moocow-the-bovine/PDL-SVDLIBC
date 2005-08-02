#!/usr/bin/perl -wd

use lib qw(./blib/lib ./blib/arch);
use PDL;
use PDL::SVDLIBC;
use PDL::CCS;
use PDL::MatrixOps;

use Benchmark qw(timethese cmpthese);

BEGIN{ $, = ' '; $eps=1e-6; }

##---------------------------------------------------------------------
## test: data
use vars qw($p $ptr $rowids $nzvals);
sub tdata {
  $p = pdl(double, [
		    [10,0,0,0,-2,0,0],
		    [3,9,0,0,0,3,1],
		    [0,7,8,7,0,0,0],
		    [3,0,8,7,5,0,1],
		    [0,8,0,9,9,13,0],
		    [0,4,0,0,2,-1,1],
		   ]);
  udata();
}

##-- gendata($n,$m,$density) : generate random data
sub gendata {
  ($n,$m,$density,$mult)=@_;
  $density = .5 if (!defined($density));
  $mult    = 1  if (!defined($mult));
  $p = random(double,$n,$m);
  $p->where($p>$density) .= 0;
  $p *= $mult;
  udata();
}

##-- gendata($n,$m,$density) : generate random data
sub gendatag {
  ($n,$m,$density,$mult)=@_;
  $density = .5 if (!defined($density));
  $mult    = 1  if (!defined($mult));
  $p = grandom(double,$n,$m);
  $p->where($p>$density) .= 0;
  $p *= $mult;
  udata();
}

##-- update cccs on changed $p
sub udata {
  ($n,$m) = $p->dims;
  ($ptr,$rowids,$nzvals) = ccsencode($p);
  ##-- HACK: allocate an extra slot in $ptr
  $ptr->reshape($ptr->nelem+1);
  $ptr->set(-1, $rowids->nelem);
}

use vars qw($ptr1 $rowids1 $nzvals1);
sub tccs {
  ($n,$m) = $p->dims;
  $nnz = $p->flat->nnz;
  svdccsencode($p,
	       ($ptr1=zeroes(long,$n)-1),
	       ($rowids1=zeroes(long,$nnz)-1),
	       ($nzvals1=zeroes(long,$nnz)-1));
}

use vars qw($iters $end $kappa $ut $ul $ssl $vt $vl $pbr $pbri);
sub tlas2 {
  #$d=4;
  $d = $p->dim(0);
  svdlas2($ptr,$rowids,$nzvals, $m,
	  ($iters=pdl(14)),
	  ($end=pdl(double,[-1e-30,1e-30])),
	  ($kappa=pdl(1e-6)),
	  ($ult = $ut = zeroes(double,$m,$d)),
	  ($sl  = $s  = zeroes(double,$d)),
	  ($vlt = $vt = zeroes(double,$n,$n)), ##-- $d-by-$n
	 );

  ##-- stretch, tranpose
  $ul  = $u  = $ult->xchg(0,1);
  $ssl = $ss = stretcher($sl);
  $vl  = $v  = $vlt->xchg(0,1);
}


use vars qw($udt $ud $sd $ssd $vdt $vd);
sub tlas2d {
  $d = $p->dim(0);
  svdlas2d($p,
	  ($iters=pdl(14)),
	  ($end=pdl(double,[-1e-30,1e-30])),
	  ($kappa=pdl(1e-6)),
	  ($udt = $ut = zeroes(double,$m,$d)),
	  ($sd  = $s  = zeroes(double,$d)),
	  ($vdt = $vt = zeroes(double,$n,$n)),
	 );

  ##-- stretch, tranpose
  $ud  = $u  = $udt->xchg(0,1);
  $vd  = $v  = $vdt->xchg(0,1);
  $ssd = $ss = stretcher($sd);
}

sub tbuiltin {
  ##-- test: compare w/ builtin svd
  ($ub,$sb,$vb) = svd($p);
  $ssb = stretcher($sb);
}


##-- check copy-ability
sub tcheck {
  my $label = shift;
  $label = '(nolabel)' if (!$label);

  ##-- hack
  $u->inplace->setnantobad->inplace->setbadtoval(0);
  $v->inplace->setnantobad->inplace->setbadtoval(0);
  $ut->inplace->setnantobad->inplace->setbadtoval(0);
  $vt->inplace->setnantobad->inplace->setbadtoval(0);

  ##-- test: copy
  $p2  = $u x $ss x $v->xchg(0,1);

  ##-- check
  print "$label : ", (all($p2->approx($p),$eps) ? "ok" : "NOT ok"), "\n";
}

sub checkall {
  tlas2;
  tcheck('las2');

  tbuiltin;
  tcheck('builtin');

  tlas2d;
  tcheck('las2d');
}

#tdata; tlas2;




##--------------------------------------------------------------------
## Restricted

use vars qw($ubr $sbr $ssbr $vbr $ur $sr $ssr $vr $pr $pri);
sub tbuiltinr {
  $dr=4 if (!defined($dr));
  $d=$dr;
  $d_1 = $d-1;

  ($ub,$sb,$vb)=($u,$s,$v)=svd($p);
  $ssb = $ss = stretcher($sb);

  ##-- restrict
  $ubr  = $ur  = $ub->slice("0:$d_1,:");
  $sbr  = $sr  = $sb->slice("0:$d_1");
  $ssbr = $ssr = stretcher($sbr);
  $vbr  = $vr  = $v->slice("0:$d_1,:");

  ##-- apply restriction
  #$pbr  = $pr  = $p x $vr;
  #$pbri = $pri = ($vr x $ssr x $ur->xchg(0,1))->xchg(0,1);
}


use vars qw($ulr $slr $sslr $vlr);
sub tlas2r {
  $dr=4 if (!defined($dr));
  $d=$dr;
  $d_1 = $d-1;

  svdlas2($ptr,$rowids,$nzvals, $m,
	  ($iters=pdl(14)),
	  ($end=pdl(double,[-1e-30,1e-30])),
	  ($kappa=pdl(1e-6)),
	  ($ultr=zeroes(double,$m,$d)),
	  ($slr=zeroes(double,$d)),
	  ($vlt=zeroes(double,$n,$n)), ##-- $n-by-$d
	 );

  ##-- stretch, tranpose
  $ulr  = $ur  = $ultr->xchg(0,1);
  $vlr  = $vr  = $vlt->slice(":,0:$d_1")->xchg(0,1);
  $sslr = $ssr = stretcher($slr);

  ##-- apply restriction
  #$plr  = $pr  = $p x $vlr;
  #$plri = $pri = ($vr x $ssr x $ur->xchg(0,1))->slice(":,:,(0)")->xchg(0,1);
}

sub tlas2dr {
  $dr=4 if (!defined($dr));
  $d=$dr;
  $d_1 = $d-1;

  svdlas2d($p,
	   ($iters=pdl(14)),
	   ($end=pdl(double,[-1e-30,1e-30])),
	   ($kappa=pdl(1e-6)),
	   ($ultr=zeroes(double,$m,$d)),
	   ($slr=zeroes(double,$d)),
	   ($vlt=zeroes(double,$n,$n)), ##-- $n-by-$d
	  );

  ##-- stretch, tranpose
  $ulr  = $ur  = $ultr->xchg(0,1);
  $vlr  = $vr  = $vlt->slice(":,0:$d_1")->xchg(0,1);
  $sslr = $ssr = stretcher($slr);

  ##-- apply restriction
  #$plr  = $pr  = $p x $vlr;
  #$plri = $pri = ($vr x $ssr x $ur->xchg(0,1))->slice(":,:,(0)")->xchg(0,1);
}

##-- checkr($label)
##   + uses current values of $dr, $d_1, $ur, $ssr, $vr
sub checkr {
  my ($label)=@_;
  $label = '(nolabel)' if (!$label);

  ##-- apply restriction
  $plr  = $pr  = $p x $vr;
  $plri = $pri = ($vr x $ssr x $ur->xchg(0,1))->slice(":,:,(0)")->xchg(0,1);

  print
    ("$label : decomp(",
     (all($plr->approx($ur x $ssr)) ? "ok" : "NOT ok"),
     "), ",
     "avg(err)=",sprintf("%.2g", abs($plri-$p)->avg), "\n",
    );
}

use vars qw($ut0 $s0 $ss0 $vt0 $utr $sr $ssr $vtr $vtrtmp $utr0 $sr0 $ssr0 $vtr0);
sub checklas2r {
  $dr = 4 if (!defined($dr));
  $drm1 = $dr-1;
  $d=4;

  tlas2;
  ($ut0,$s0,$ss0,$vt0)=($ut,$s,$ss,$vt);
  $utr0 = $ut->slice(":,0:$drm1");
  $sr0  = $s->slice("0:$drm1");
  $ssr0 = stretcher($sr0);
  $vtr0 = $vt->slice("0:$drm1,:");

  tlas2r;
  ($utr,$sr,$ssr,$vtrtmp)=($ut,$s,$ss,$vt);
  $vtr  = $vtrtmp->slice("0:$drm1,:");

  ##-- copy 'em
  #$p0c = $utr0->xchg(0,1) x $ssr0 x $vtr0;
}


##---------------------------------------------------------------------
## DUMMY
##---------------------------------------------------------------------
foreach $i (0..100) {
  print "--dummy($i)--\n";
}

