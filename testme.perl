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
use vars qw($p $ptr $rowids $nzvals $a);
sub tdata {
  $p = $a = pdl(double, [
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
  ##-- ccs encode
  ($ptr,$rowids,$nzvals) = ccsencode($p);
  #$ptr = $ptr->convert(longlong);
  #$rowids = $rowids->convert(longlong);
  ##-- HACK: allocate an extra slot in $ptr
  $ptr->reshape($ptr->nelem+1);
  $ptr->set(-1, $rowids->nelem);
}

use vars qw($ptr1 $rowids1 $nzvals1);
sub tccs {
  ($n,$m) = $p->dims;
  $nnz = $p->flat->nnz;
  _svdccsencode($p,
		($ptr1=zeroes(long,$n+1)-1),
		($rowids1=zeroes(long,$nnz)-1),
		($nzvals1=zeroes(long,$nnz)-1));
}

use vars qw($iters $end $kappa $ut $ul $ssl $vt $vl $pbr $pbri);
sub tlas2 {
  #$d=4;
  tdata() if (!defined($p));

  $d = $p->dim(0);
  svdlas2($ptr,$rowids,$nzvals, $m,
	  ($iters=pdl(14)),
	  ($end=pdl(double,[-1e-30,1e-30])),
	  ($kappa=pdl(1e-6)),
	  ($u2t = $u2t = zeroes(double,$m,$d)),
	  ($s2  = $s2  = zeroes(double,$d)),
	  ($v2t = $v2t = zeroes(double,$n,$d)), ##-- $d-by-$n
	 );

  ##-- stretch, tranpose
  use vars qw($s2s);
  $u2  = $u  = $u2t->xchg(0,1);
  $s2s = $ss = stretcher($s2);
  $v2  = $v  = $v2t->xchg(0,1);
}
#tlas2();


use vars qw($udt $ud $sd $ssd $vdt $vd);
sub tlas2d { #-- ok
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
  my ($label,$u,$s,$v) = @_;
  $label = '(nolabel)' if (!$label);

  ##-- hack
  $u->inplace->setnantobad->inplace->setbadtoval(0);
  $s->inplace->setnantobad->inplace->setbadtoval(0);
  $v->inplace->setnantobad->inplace->setbadtoval(0);

  ##-- test: copy
  $p2  = $u x stretcher($s) x $v->xchg(0,1);

  ##-- check
  print "$label : ", (all($p2->approx($p),$eps) ? "ok" : "NOT ok"), "\n";
}

sub checkall {
  tdata() if (!defined($p));

  tbuiltin;
  tcheck('builtin', $ub,$sb,$vb);

  tlas2d;
  tcheck('las2d', $ud, $sd, $vd);

  tlas2; ##-- strangeness with SVDLIBC 'long' vs. PDL 'int' or 'longlong' (now using: int)
  tcheck('las2', $u2, $s2, $v2);
}
#checkall;




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
}


use vars qw($ulr $slr $sslr $vlr);
sub tlas2dr {
  tdata() if (!defined($p));
  $dr=4 if (!defined($dr));
  $d=$dr;
  $d_1 = $d-1;

  svdlas2d($p,
	   ($iters=pdl(14)),
	   ($end=pdl(double,[-1e-30,1e-30])),
	   ($kappa=pdl(1e-6)),
	   ($udtr=zeroes(double,$m,$d)),
	   ($sdr=zeroes(double,$d)),
	   ($vdtr=zeroes(double,$n,$d)), ##-- $n-by-$d
	  );

  ##-- stretch, tranpose
  $udr  = $ur  = $udtr->xchg(0,1);
  $vdr  = $vr  = $vdtr->xchg(0,1);
  $ssdr = $ssr = stretcher($sdr);

  ##-- apply restriction
  #$plr  = $pr  = $p x $vlr;
  #$plri = $pri = ($vr x $ssr x $ur->xchg(0,1))->slice(":,:,(0)")->xchg(0,1);
}

sub tlas2r {
  tdata if (!defined($p));
  $dr=4 if (!defined($dr));
  $d=$dr;
  $d_1 = $d-1;

  svdlas2($ptr,$rowids,$nzvals, $m,
	  ($iters=pdl(14)),
	  ($end=pdl(double,[-1e-30,1e-30])),
	  ($kappa=pdl(1e-6)),
	  ($u2rt=zeroes(double,$m,$d)),
	  ($s2r=zeroes(double,$d)),
	  ($v2rt=zeroes(double,$n,$d)), ##-- $n-by-$d
	 );

  ##-- stretch, tranpose
  $u2r  = $ur  = $u2rt->xchg(0,1);
  $v2r  = $vr  = $v2rt->xchg(0,1);
  $ss2r = $ssr = stretcher($s2r);

  ##-- apply restriction
  #$plr  = $pr  = $p x $vlr;
  #$plri = $pri = ($vr x $ssr x $ur->xchg(0,1))->slice(":,:,(0)")->xchg(0,1);
}


##-- tcheckr($label, $ur,$sr,$vr)
##   + uses current values of $dr, $d_1, $ur, $ssr, $vr  : ?!?!?!
use vars qw($ss2r $ssdr);
sub tcheckr {
  my ($label,$ur,$sr,$vr)=@_;
  $label = '(nolabel)' if (!$label);

  ##-- bad-value hack
  $ur->inplace->setnantobad->inplace->setbadtoval(0);
  $sr->inplace->setnantobad->inplace->setbadtoval(0);
  $vr->inplace->setnantobad->inplace->setbadtoval(0);

  ##-- apply restriction
  $plr  = $pr  = $p x $vr;
  $plri = $pri = ($vr x stretcher($sr) x $ur->xchg(0,1))->slice(":,:,(0)")->xchg(0,1);

  print
    ("$label : decomp(",
     (all($plr->approx($ur x stretcher($sr))) ? "ok" : "NOT ok"),
     "), ",
     "avg(err)=",sprintf("%.2g", abs($plri-$p)->avg), "\n",
    );
}


sub checkallr {
  tdata() if (!defined($p));

  tbuiltinr;
  tcheckr('builtin', $ubr,$sbr,$vbr);

  tlas2dr;
  tcheckr('las2d', $udr, $sdr, $vdr);

  tlas2r;
  tcheckr('las2', $u2r, $s2r, $v2r);
}
checkallr;


##---------------------------------------------------------------------
## DUMMY
##---------------------------------------------------------------------
foreach $i (0..3) {
  print "--dummy($i)--\n";
}

