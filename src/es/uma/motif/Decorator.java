package es.uma.motif;

class Decorator implements Comparable<Object> {
    int init, end;
    String classname;
    boolean[] mask;

    Decorator(int start, int end, String classname, boolean[] mask) {
        this.init = start;
        this.end = end;
        this.classname = classname;
        this.mask = mask;
    }

    @Override
    public int compareTo(Object o) {
        Decorator target = (Decorator) o;
        return (init == target.init) ? 0 : ((init < target.init)? -1 : 1);
    }
}
