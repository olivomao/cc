#import tensorflow as tf
import numpy as np
np.random.seed(0)
import random
random.seed(0)
import editdistance
import pdb
import matplotlib #draw_hist issue
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from cgk import proj_hamming_dist

simType = ['sub', 'ins', 'del', 'match']
simTypeSize = 4

'''
items [[val1_list, val1_description],
       [val2_list, val2_description],
       ...]

output to output_fig

draw hist plot where per curve corresponds
to the histogram of val_list, and is labeled by
 val_description

'''
def draw_histogram(items,
                   output_fig):

    fig, ax = plt.subplots()

    for itm in items:
        v_list = itm[0]
        v_description = itm[1]
        ax.hist(v_list, 50, 
                histtype='step', 
                normed=True,
                label=v_description)

    ax.legend()
    plt.tight_layout()
    plt.savefig(output_fig)
    print('draw_histogram %s drawn'%output_fig)
    return

'''
a distance class that calculates dist b/w
two strings a, b

the dist calc includes:
- edit distance
- misc (e.g. in cgk.py) ==> to be incorporated
'''
class Distance(object):

    def __init__(self):
        print('Init Distance')

        return

    def edit_dist(self,a,b):
        #return float(Levenshtein.distance(a,b))/(0.5*(len(a)+len(b)))
        #return float(editdistance.eval(a,b))/(0.5*(len(a)+len(b)))
        #return float(editdistance.eval(a,b))/(len(a)+len(b))
        return editdistance.eval(a,b)

'''
an indel channel
'''
class Channel(object):

    def __init__(self,
                 ins_rate=0.0,
                 del_rate=0.0,
                 sub_rate=0.0):

        self.ins_rate = ins_rate
        self.del_rate = del_rate
        self.sub_rate = sub_rate

        print('Init Channel (i=%f, d=%f, s=%f)'%\
              (self.ins_rate,
               self.del_rate,
               self.sub_rate))
        return

    '''
    input s goes through the channel and get
    output t
    '''
    def output(self, s):
        t = self.sim_seq_binary(s,
                                sub=self.sub_rate,
                                dele=self.del_rate,
                                ins=self.ins_rate)
        return t 

    '''
    remove a block of length L from s (L<|s|)
    the block can start at any pos (<=Ls-L) of s
    '''
    def rand_remove_block(self, s, L):
        assert len(s)>L
        i = random.randint(0, len(s)-L) #[0, Ls-L]
        t = s[0:i]+s[i+L:len(s)]

        print('rand remove block: s=%s i=%d t=%s'%(s, i, t))

        return t

    '''
    add a block of length L into s
    the block can start at any pos of s
    '''
    def rand_add_block(self, s, L):
        i = random.randint(0, len(s)-1)
        d = gen_s(L)
        t = s[0:i]+d+s[i:len(s)]

        print('rand add block: s=%s i=%d d=%s t=%s'%(s, i, d, t))
        #pdb.set_trace() 
        return t

    def rand_flip_block(self, s, L):
        def flip(st):
            return ''.join('1' if x == '0' else '0' for x in st)
        assert len(s)>L
        i = random.randint(0, len(s)-L) #[0, Ls-L]
        t = s[0:i]+flip(s[i:i+L])+s[i+L:len(s)] 
        print('rand flip block: s=%s i=%d t=%s'%(s, i, t))
        return t

    '''
    depending on tp ('add', 'remove', 'flip')
    rand add/remove/flip a block of len L starting from any
    possible loc of s
    '''
    def rand_block(self, s, L, tp):
        if tp == 'add':
            t = self.rand_add_block(s, L)
        elif tp == 'remove':
            t = self.rand_remove_block(s, L)
        elif tp == 'flip':
            t = self.rand_flip_block(s, L)
        else:
            print('rand_block unknown type: %s'%tp)
            pdb.set_trace()
            t = None
        return t
    

    '''
    simulate a binary seq through a long read channel

    input:
    seq - input binary seq
    sub - sub rate
    ins - ins rate
    dele - dele rate
    (need: 0 <= sub + ins + del < 1)
    '''
    def sim_seq_binary(self, seq, sub, dele, ins):

        if(sub+ins+dele>=1):
            print >> sys.stderr, "Total sub+ins+del error can not exceed 1!"
            sys.exit(-1)

        profile = [sub, sub+ins, sub+ins+dele, 1.]

        return self.sim_seq_binary_profile(seq, profile)
        

    def sim_seq_binary_profile(self, seq, profile):

        nucl = set(['0', '1'])
        sim = ''
        for i, s in enumerate(seq):
            while True:
                tp = self.throwdice(profile)
                if tp=='match': 
                    sim += s
                    break
                elif tp=='ins':
                    # insertion occurred, with 1/4 chance it'll match
                    choice = random.sample(nucl,1)[0]
                    sim += choice
                elif tp=='sub': # anything but the right one
                    choice = random.sample(nucl.difference([s]),1)[0]
                    sim += choice
                    break
                elif tp=='del': # skip over this
                    break
                else: 
                    print >> sys.stderr, "Invalid type %s"%(tp)
                    sys.exit(-1)


        return sim

    def throwdice(self, profile):
        dice = random.random()
        for i in xrange(simTypeSize):
            if dice < profile[i]:
                return simType[i]

'''
generate random bin string s of length L
'''
def gen_s(L):
    letters = "01"
    s = "".join(random.choice(letters) for i in range(L))
    return s

'''
return R_list to be used by proj_hamming_dist
'''
def init_R_list(L, N_iter=1):
    L_R = 3*L*2

    R_list = []
    for i in range(N_iter):
        R = np.random.randint(0,2,L_R)
        R_list.append(R)

    return R_list

def main():

    print('main')
    ######## configs
    N = 2000
    L = 100
    ch = Channel(ins_rate=0.0,
                 del_rate=0.0, #0.0,
                 sub_rate=0.3) #0.3)
    t2_type = 'flip' #remove, add, flip
    dist = Distance()

    #fix a rand matrix
    R_list = None #init_R_list(L, N_iter=1) 

    fix_sample = True #False # True
    ######## configs end
 
    '''
    #tensorboard utilization (does not seem to work now)
    de_s_t1 = tf.get_variable(name='de_s_t1', shape=[N], trainable=False)
    tf.summary.histogram('de_s_t1', de_s_t1)
    merged = tf.summary.merge_all()
    writer = tf.summary.FileWriter('/data1/shunfu1/clustering_cgk_vs_edit/log/')
    '''

    d_s_t1_acc = []
    d_s_t2_acc = []
    d_cgk_s_t1_acc = []
    d_cgk_s_t2_acc = []

    if fix_sample == True:
        fix_s = gen_s(L)
        fix_t1 = ch.output(fix_s)
        fix_t2 = ch.rand_block(fix_s,
                               dist.edit_dist(fix_s, fix_t1),
                               t2_type)

    for i in range(N):
        print(i)

        if fix_sample == True:
            s = fix_s
        else:
            s = gen_s(L)
        print('s='+s)

        if fix_sample == True:
            t1 = fix_t1
        else:
            t1 = ch.output(s)
        print('t1='+t1)

        d_s_t1 = dist.edit_dist(s,t1)
        d_s_t1_acc.append(d_s_t1)
        print('edit_dist(s,t1)=%d'%d_s_t1)

        if fix_sample == True:
            t2 = fix_t2
        else:
            t2 = ch.rand_block(s, d_s_t1, t2_type)
        print('t2='+t2)
        d_s_t2 = dist.edit_dist(s,t2)
        d_s_t2_acc.append(d_s_t2)
        print('edit_dist(s,t2)=%d'%d_s_t2)

        d_cgk_s_t1 = proj_hamming_dist(s, t1, 0, False, R_list)
        d_cgk_s_t1_acc.append(d_cgk_s_t1)

        d_cgk_s_t2 = proj_hamming_dist(s, t2, 0, False, R_list)
        d_cgk_s_t2_acc.append(d_cgk_s_t2)

        print('')

    draw_histogram([[d_s_t1_acc, 'de_s_t1'],
                    [d_s_t2_acc, 'de_s_t2']], 
                   'fig.png') 

    draw_histogram([[d_cgk_s_t1_acc, 'de_cgk_s_t1'],
                    [d_cgk_s_t2_acc, 'de_cgk_s_t2']], 
                   'fig_cgk.png') 
    '''
    # tensorbloard utilization (does not seem to work)
    summ, d_ = sess.run([merged, 
                        de_s_t1],
                        feed_dict={de_s_t1: np.asarray(d_acc)})
    print('d_=%s'%str(d_))
    print('')
    writer.add_summary(summ)#all samples considered at time 0
    '''

    return


if __name__ == "__main__":

    main()


